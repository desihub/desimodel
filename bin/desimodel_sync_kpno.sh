#!/bin/bash

# Exit on error
set -e

# This script is intended to be run from a cronjob (see 
# desimodel/etc/desimodel_sync_kpno_cron.sh).  You should not run
# this script interactively.  It expects to be run inside an
# already-configured software environment.

fpsync=$(which desi_sync_focalplane)
echo "Using focalplane sync script:  ${fpsync}"

# Find the newest calibration file
caldir="/data/focalplane/calibration/"
calfile=$(ls ${caldir} | egrep '[0-9]{8}T[0-9]{6}.*' | sort | tail -n 1)
calpath="${caldir}/${calfile}"

# Run it.
eval ${fpsync} --calib_file ${calpath} --commit
