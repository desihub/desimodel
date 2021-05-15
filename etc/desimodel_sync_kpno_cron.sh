#!/bin/bash

# Exit on error
set -e

# Path to the sync directory, containing the desimodel data trunk
# checkout and logs.
syncdir=/n/home/datasystems/desimodel_sync
logdir="${syncdir}/logs"
svntrunk="${syncdir}/desimodel_trunk"

logdate="$(date -u --iso-8601=seconds)"
logname="sync_${logdate}"
logfile="${logdir}/${logname}"

echo "Running at ${logdate}" > "${logfile}"

# Get the latest desiconda install
desiconda=$(ls -d /software/datasystems/desiconda/20* | sort | tail -n 1)

# Get the latest KPNO version of desimodules
desimodules=$(ls ${desiconda}/modulefiles/desimodules/*-kpno | sort | tail -n 1 | xargs basename)

echo "Using latest desiconda:  ${desiconda}" >> "${logfile}"
echo "Using latest KPNO version of desimodules:  ${desimodules}" >> "${logfile}"

# Set up environment

export DESI_PRODUCT_ROOT="${desiconda}"
export DESI_ROOT=/data/datasystems
export DESI_TARGET=${DESI_ROOT}/target
export DESI_SURVEYOPS=${DESI_ROOT}/survey/ops/surveyops/trunk

module use ${DESI_PRODUCT_ROOT}/modulefiles
module load desiconda
module load desimodules/${desimodules}

echo "Using desimodel data svn trunk at ${svntrunk}" >> "${logfile}"
export DESIMODEL="${svntrunk}"

# Find our sync script and run it.
runsync=$(which desimodel_sync_kpno.sh)
echo "Calling sync script:  ${runsync}" >> "${logfile}"

eval ${runsync} >> "${logfile}"

echo "Sync script finished at $(date -u --iso-8601=seconds)" >> "${logfile}"
