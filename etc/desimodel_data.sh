#!/bin/bash -x
set -e
svn export https://desi.lbl.gov/svn/code/desimodel/${DESIMODEL_VERSION}/data
[[ -z "${DESIMODEL}" ]] && export DESIMODEL=$(pwd)
