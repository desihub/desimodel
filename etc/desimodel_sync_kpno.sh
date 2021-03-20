#!/bin/bash

# Exit on error
set -e

# Location of the minimal virtualenv with just desiutil / desimodel
# and dependencies.
venv_sync=/software/datasystems/users/kisner/desimodel_sync

# Location of logs
sync_logs=/software/datasystems/users/kisner/logs

# Location of desimodel trunk and "for_sims" branches
desimodel_trunk=/software/datasystems/users/kisner/svn/desimodel_trunk
desimodel_sim=/software/datasystems/users/kisner/svn/desimodel_for_sims

# Load EUPS
source /software/products/eups-2.1.4/bin/setups.sh

# Load default python3
setup python

# NOTE: create a minimal python environment with (for example)
# python3 -m venv ${venv_sync}

# Activate our minimal python environment
source ${venv_sync}/bin/activate

# NOTE: on first use, install our packages:
# python3 -m pip install --upgrade pip
# python3 -m pip install wheel
# python3 -m pip install pyyaml scipy matplotlib numpy astropy configobj xlrd requests
# python3 -m pip install 'git+https://github.com/desihub/desiutil.git@master#egg=desiutil'
# python3 -m pip install 'git+https://github.com/desihub/desimodel.git@sync_focalplane#egg=desimodel'

# Find the most recent DB dump file
caldir="/data/focalplane/calibration/"
calfile=$(ls ${caldir} | egrep '[0-9]{8}T[0-9]{6}.*' | sort | tail -n 1)
calpath="${caldir}/${calfile}"

# Log root
log="${sync_logs}/sync_$(date "+%Y%m%d-%H%M%S")"

# Run the sync for the main (trunk) data location
DESIMODEL="${desimodel_trunk}" desi_sync_focalplane \
    --calib_file ${calpath} \
    --commit --test 2>&1 > "${log}_trunk"

# Run the sync for the "simulation" branch, which has all transient states
# (restricted reach) cleared.
DESIMODEL="${desimodel_sim}" desi_sync_focalplane \
    --calib_file ${calpath} \
    --simulate_good \
    --commit --test 2>&1 > "${log}_sim"
