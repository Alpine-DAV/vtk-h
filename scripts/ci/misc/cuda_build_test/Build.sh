#!/bin/bash
set -e
date
#################################
# setup build env
#################################
# output env
cat etc/*rel*
env
whoami

export BRANCH=$1
echo "BRANCH: $BRANCH"
#
# check out branch we want to build
#
rm -rf vtk-h
git clone --recursive --depth=1 https://github.com/Alpine-DAV/vtk-h.git --branch=${BRANCH}
cd vtk-h

# use updated config
cp /config.yaml scripts/uberenv/spack_configs/ci/ubuntu_16_cuda_10.1_devel/
echo "Spack Config:"
cat scripts/uberenv/spack_configs/ci/ubuntu_16_cuda_10.1_devel/config.yaml

#################################
# run uber to build tpls
#################################
# echo system python details
which python
python --version
# setup spack spec (must be static)
export SPACK_SPEC="%gcc+mpi+cuda~openmp~shared"
export SPACK_SPEC="${SPACK_SPEC} ^vtk-m+cuda~shared"
export SPACK_SPEC="${SPACK_SPEC} ^cmake~openssl~ncurses"
echo $SPACK_SPEC
# run uber to build tpls
date
python scripts/uberenv/uberenv.py --pull --spec "${SPACK_SPEC}" --spack-config-dir=scripts/uberenv/spack_configs/ci/ubuntu_16_cuda_10.1_devel/
date
#################################
# configure
#################################
# setup compiler env vars
export CC=gcc
export CXX=g++
export FC=gfortran
${CC} --version
# capture current path
export ROOT_DIR=`pwd`
# find spack generated host config file
export HOST_CONFIG=`ls ${ROOT_DIR}/uberenv_libs/*.cmake`
echo "HOST_CONFIG: $HOST_CONFIG"
# find spack installed cmake
export CMAKE_BIN_DIR=`ls -d ${ROOT_DIR}/uberenv_libs/spack/opt/spack/*/*/cmake*/bin`
export PATH=${CMAKE_BIN_DIR}:$PATH
echo $PATH
which cmake
cmake --version
# prepare build dir
mkdir build
cd build
# setup cmake options
export CMAKE_OPTS="-DCMAKE_BUILD_TYPE=Release"
export CMAKE_OPTS="${CMAKE_OPTS} -DBUILD_SHARED_LIBS=OFF"
export CMAKE_OPTS="${CMAKE_OPTS} -DCMAKE_INSTALL_PREFIX=../install"
# configure
cmake ${CMAKE_OPTS} -C ${HOST_CONFIG} ../src
#################################
# build
#################################
# build
date
time make -j 20 VERBOSE=1
date
#################################
# install
#################################
make install
#################################
# check install
#################################
cd ..
ls install
date