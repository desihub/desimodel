#!/bin/bash

# IMPORTANT:  If you update this script, remember to copy it to
# ~datasystems/desimodel_sync/ at KPNO.  This copy is needed so that
# cron can find it in a fixed location independent of particular
# versions of the software stack (this script loads the software
# stack).

# Ensure that we are running this as the datasystems user
if [ $(whoami) != "datasystems" ]; then
    echo "You should only run this script at KPNO as the datasystems user"
    exit 1
fi

# Path to the sync directory, containing the desimodel data trunk
# checkout and logs.
syncdir=/n/home/datasystems/desimodel_sync
logdir="${syncdir}/logs"
svntrunk="${syncdir}/desimodel_trunk"

logdate="$(date -u --iso-8601=seconds)"
logname="sync_${logdate}"
logfile="${logdir}/${logname}"

# Default sync status until confirmed otherwise
sync_status="FAILED"

echo "Running at ${logdate}" > "${logfile}"

# Use the latest default desiconda version
desiconda=/software/datasystems/desiconda/kpno-20250320-2.2.1.dev

# Get the latest stable version of desimodules
desimodules=$(ls -d ${desiconda}/modulefiles/desimodules/2* | sort -V | tail -n 1 | xargs basename)

echo "Using default desiconda:  ${desiconda}" >> "${logfile}"
echo "Using latest stable version of desimodules:  ${desimodules}" >> "${logfile}"

# Set up environment
export DESI_PRODUCT_ROOT="${desiconda}"
export DESI_ROOT=/data/datasystems
export DESI_TARGET=${DESI_ROOT}/target
export DESI_SURVEYOPS=${DESI_ROOT}/survey/ops/surveyops/trunk
export DESIMODEL_CENTRAL_REPO=${DESI_ROOT}/survey/ops/desimodel/trunk

module use ${DESI_PRODUCT_ROOT}/modulefiles
module load desiconda
module load desimodules/${desimodules}

# 2025-11-19: this hack is needed by cron to set up the paths. Investigate...
module swap -f desiutil/3.5.0
module swap -f desitree/0.7.0
module swap -f specter/0.10.1
module swap -f gpu_specter/0.2.1
module swap -f desimodel/0.19.3
module swap -f desitarget/2.9.0
module swap -f specsim/v0.17
module swap -f desispec/0.69.0
module swap -f desisim/0.38.1
module swap -f fiberassign/5.7.2
module swap -f desisurvey/0.20.0
module swap -f surveysim/0.12.6
module swap -f redrock/0.20.4
module swap -f redrock-templates/0.9.1
module swap -f desimeter/0.7.1
module swap -f simqso/v1.3.0
module swap -f speclite/v0.20

echo "Using desimodel data svn trunk at ${svntrunk}" >> "${logfile}"
export DESIMODEL="${svntrunk}"

# Find the sync script
fpsync=$(which desi_sync_focalplane)
echo "Using focalplane sync script:  ${fpsync}" >> "${logfile}"

# Find the newest calibration file
caldir="/data/focalplane/calibration/"
calfile=$(ls ${caldir} | egrep '[0-9]{8}T[0-9]{6}.*' | sort | tail -n 1)
calpath="${caldir}${calfile}"
echo "Found newest calibration file:  ${calpath}" >> "${logfile}"

# Make sure that any locally modified files are removed
echo "Ensuring clean svn tree at ${DESIMODEL}" >> "${logfile}"
svn revert -R "${svntrunk}/data" >> "${logfile}"
svn up "${svntrunk}/data" >> "${logfile}"

# Run it, without committing result.
eval ${fpsync} --calib_file ${calpath} >> "${logfile}" 2>&1
syncstatus=$?

# Check for clean execution and no errors in the output log.
nerr=`grep -c ERROR ${logfile}`

if [ ${syncstatus} -ne 0 ] || [ ${nerr} -ne 0 ]; then
    echo "Focalplane sync failed" >> "${logfile}"
else
    echo "Focalplane sync completed" >> "${logfile}"
    # Make sure we can load the resulting hardware model
    echo "Try loading focalplane model at ${DESIMODEL}... " >> "${logfile}"
    PYTHON_CODE=$(cat <<END
from fiberassign.hardware import load_hardware
hw = load_hardware()
END
)
    result=$(python3 -c "$PYTHON_CODE")
    if [ $? -eq 0 ]; then
        sync_status="SUCCESSFUL"
        echo "SUCCESS" >> "${logfile}"
        # Now commit result
        mesg="Appending DB sync ${calfile} to current focalplane model"
        svn commit -m "${mesg}" "${svntrunk}/data" >> "${logfile}"
        svn up "${svntrunk}/data" >> "${logfile}"
        echo "Updating $DESIMODEL_CENTRAL_REPO." >> "${logfile}"
        svn up "${DESIMODEL_CENTRAL_REPO}" >> "${logfile}"
    else
        echo "FAIL" >> "${logfile}"
        echo "Refusing to commit broken update.  Restore the *.previous files." >> "${logfile}"
    fi
fi

# Send notifications. Mail out result.
slack_email=${DESI_SLACK_MAIL_DESIMODEL_SYNC}

if [ "x${slack_email}" = "x" ]; then
    echo "Environment variable DESI_SLACK_MAIL_DESIMODEL_SYNC not set- skipping notifications" >> "${logfile}"
else
    cat ${logfile} | mailx -s "Focalplane sync ${sync_status}: ${logdate}" ${slack_email}
    echo "Sent log to slack" >> ${logfile}
fi

echo "Sync script finished at $(date -u --iso-8601=seconds)" >> "${logfile}"

# Guard against accidental deletion
chmod a-w "${logfile}"
