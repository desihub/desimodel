#!/bin/bash -x
[[ ${DEBUG} == True ]] && set -x
export DESIMODEL=${HOME}/desimodel/${DESIMODEL_VERSION}
mkdir -p ${DESIMODEL}
cd ${DESIMODEL}
svn export https://desi.lbl.gov/svn/code/desimodel/${DESIMODEL_VERSION}/data
[[ ${DEBUG} == True ]] && set +x
