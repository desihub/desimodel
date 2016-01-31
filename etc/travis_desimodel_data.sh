#!/bin/bash -x
set -e
export DESIMODEL=${HOME}/desimodel/${DESIMODEL_VERSION}
mkdir -p ${DESIMODEL}
# Do this in a subshell, so the directory change doesn't affect the main script.
(cd ${DESIMODEL} && svn export https://desi.lbl.gov/svn/code/desimodel/${DESIMODEL_VERSION}/data)
