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

echo "Running at ${logdate}" > "${logfile}"

# Use the latest default desiconda version
desiconda=/software/datasystems/desiconda/default

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
module swap -f desimodel/0.17.0

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
if [ $? -ne 0 ]; then
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

# Send notifications.

# Get our webhook address from the environment
slack_web_hook=${DESI_SLACKBOT_DESIMODEL_SYNC}

if [ "x${slack_web_hook}" = "x" ]; then
    echo "Environment variable DESI_SLACKBOT_DESIMODEL_SYNC not set- skipping notifications" >> "${logfile}"
else
    # Create the JSON payload.
    slackjson="${logfile}_slack.json"
    headtail=15
    echo -e "{\"text\":\"Focalplane DB sync (log at \`${logfile}\`):\n\`\`\`$(head -n ${headtail} ${logfile} | sed -e "s|'|\\\'|g")\`\`\`\n(Snip)\n\`\`\`$(tail -n ${headtail} ${logfile} | sed -e "s|'|\\\'|g")\`\`\`\"}" > "${slackjson}"
    # Post it.
    slackerror=$(curl -X POST -H 'Content-type: application/json' --data "$(cat ${slackjson})" ${slack_web_hook})
    echo "Slack API post  ${slackerror}" >> "${logfile}"
fi

echo "Sync script finished at $(date -u --iso-8601=seconds)" >> "${logfile}"

# Guard against accidental deletion
chmod a-w "${logfile}"
