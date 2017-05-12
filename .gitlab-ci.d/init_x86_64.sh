#!/bin/bash

#Determine which OS you are using
if [ "$(uname)" == "Darwin" ]; then
    if [ $(sw_vers -productVersion | awk -F '.' '{print $1 "." $2}') == "10.12" ]; then
        OS=mac1012
        COMPILER_TYPE=clang
        COMPILER_VERSION=clang80
    else
        echo "Bootstrap only works on macOS Sierra (10.12)"
    fi
elif [ "$(uname)" == "Linux" ]; then
    if [ "$( cat /etc/*-release | grep Scientific )" ]; then
        OS=slc6
    elif [ "$( cat /etc/*-release | grep CentOS )" ]; then
        OS=centos7
    fi
else
    echo "UNKNOWN OS"
    exit 1
fi

#Determine is you have CVMFS installed
if [ ! -d "/cvmfs" ]; then
    echo "No CVMFS detected, please install it."
    exit 1
fi

if [ ! -d "/cvmfs/clicdp.cern.ch" ]; then
    echo "No clicdp CVMFS repository detected, please add it."
    exit 1
fi


#Determine which compiler to use
if [ -z ${COMPILER_TYPE} ]; then
    COMPILER_TYPE="gcc"
fi
if [ ${COMPILER_TYPE} == "gcc" ]; then
    COMPILER_VERSION="gcc7"
fi
if [ ${COMPILER_TYPE} == "llvm" ]; then
    COMPILER_VERSION="llvm40"
fi


#Choose build type
if [ -z ${BUILD_TYPE} ]; then
    BUILD_TYPE=opt
fi


# General variables
CLICREPO=/cvmfs/clicdp.cern.ch
BUILD_FLAVOUR=x86_64-${OS}-${COMPILER_VERSION}-${BUILD_TYPE}

#--------------------------------------------------------------------------------
#     Compiler
#--------------------------------------------------------------------------------

if [ ${COMPILER_TYPE} == "gcc" ]; then
    source ${CLICREPO}/compilers/gcc/7.1.0/x86_64-${OS}/setup.sh
fi
if [ ${COMPILER_TYPE} == "llvm" ]; then
    source ${CLICREPO}/compilers/llvm/4.0.0/x86_64-${OS}/setup.sh
fi

#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------

export CMAKE_HOME=${CLICREPO}/software/CMake/3.8.1/${BUILD_FLAVOUR}
export PATH=${CMAKE_HOME}/bin:$PATH

#--------------------------------------------------------------------------------
#     ROOT
#--------------------------------------------------------------------------------

export ROOTSYS=${CLICREPO}/software/ROOT/6.08.06/${BUILD_FLAVOUR}
export PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"

#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------

export XercesC_HOME=${CLICREPO}/${BUILD_FLAVOUR}
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"

#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------

export G4INSTALL=${CLICREPO}/software/Geant4/10.03.p01/${BUILD_FLAVOUR}
export G4LIB=$G4INSTALL/lib64/Geant4-10.3.1/
export G4ENV_INIT="${G4INSTALL}/bin/geant4.sh"
export G4SYSTEM="Linux-g++"

#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------

export BOOST_ROOT=${CLICREPO}/software/Boost/1.64.0/${BUILD_FLAVOUR}
export LD_LIBRARY_PATH="${BOOST_ROOT}/lib:$LD_LIBRARY_PATH"

#--------------------------------------------------------------------------------
#     Ninja
#--------------------------------------------------------------------------------

export Ninja_HOME=${CLICREPO}/software/Ninja/1.7.2/${BUILD_FLAVOUR}
export PATH="$Ninja_HOME:$PATH"

#--------------------------------------------------------------------------------
#     Eigen
#--------------------------------------------------------------------------------

export Eigen_HOME=${CLICREPO}/software/Eigen/3.3.0/x86_64-centos7-gcc62-opt/

