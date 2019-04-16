#!/usr/bin/env bash
PKG_NAME=mantid_total_scattering
USER=marshallmcdonnell

OS=${TRAVIS_OS_NAME}-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
conda build .
export BUILD=$(ls ${CONDA_BLD_PATH}/${OS}/${PKG_NAME}* | sed -n 's/.*${PKG_NAME}-\(.*\)-\(.*\)\.tar.bz2/\1-\2/p'
anaconda -t ${CONDA_UPLOAD_TOKEN} upload -u ${USER} -l nightly ${CONDA_BLD_PATH}/${OS}/${PKG_NAME}-${BUILD}.tar.bz2 --force

