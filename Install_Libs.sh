#!/bin/sh
#set -e

rcFile="${HOME}/.bashrc"

BASEFOLDER=$(pwd)

LIBRARY_DIR=$BASEFOLDER/Libraries
if ! grep -q "^export LIBRARY_DIR=${LIBRARY_DIR}" "$rcFile"; then
    echo "export LIBRARY_DIR=${LIBRARY_DIR}" >> "$rcFile"
    echo "[inserted] export LIBRARY_DIR=${LIBRARY_DIR}"
fi

## pull ConcreteCorrosion
if [ ! -d "ConcreteCorrosion" ]; then
    git clone https://github.com/T-Hageman/ConcreteCorrosion.git
fi

#check compilers
C_COMPILER=$(which gcc)
CPP_COMPILER=$(which g++)
F_COMPILER=$(which gfortran)
if [ -z "$C_COMPILER" ] || [ -z "$CPP_COMPILER" ] || [ -z "$F_COMPILER" ]; then
    sudo apt-get -y update
    sudo apt-get -y upgrade
    sudo apt-get -y install firefox
    sudo apt-get -y install gedit
    sudo apt-get -y install build-essential
    sudo apt-get -y install gcc
    sudo apt-get -y install g++
    sudo apt-get -y install gfortran
    sudo apt-get -y install gdb
    sudo apt-get -y install cmake-curses-gui
    sudo apt-get -y install mesa-common-dev
    sudo apt-get -y install mesa-utils
    sudo apt-get -y install freeglut3-dev
    sudo apt-get -y install ninja-build
    sudo apt-get -y install pkgconf
fi

# check cmake
C_MAKE=$(which cmake)
if [ -z "$C_MAKE" ]; then
    sudo apt-get -y install cmake
fi

## Install libraries
if [ ! -d $LIBRARY_DIR ]; then
    mkdir $LIBRARY_DIR
fi

# intel oneAPI
if [ -z "$MKLROOT" ]; then
    MKL_INSTALL_FOLDER="${BASEFOLDER}/Libraries/Intel"
    sudo apt-get install xdg-utils
    sudo apt-get install libnss3

    wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh
    sudo sh ./l_BaseKit_p_2023.2.0.49397.sh -a -s --eula accept --install-dir $MKL_INSTALL_FOLDER
    rm l_BaseKit_p_2023.2.0.49397.sh

    chmod +x ${MKL_INSTALL_FOLDER}/setvars.sh
    source ${MKL_INSTALL_FOLDER}/setvars.sh intel64

    # add to bashrc
    if ! grep -q "^source ${MKL_INSTALL_FOLDER}/setvars.sh intel64" "$rcFile"; then
        echo "source ${MKL_INSTALL_FOLDER}/setvars.sh intel64" >> "$rcFile"
        echo "[inserted] source ${MKL_INSTALL_FOLDER}/setvars.sh intel64"
    fi
fi

# PETSC
if [ -z "$PETSC_DIR" ] || [ -z "$PETSC_ARCH" ]; then
    PETSC_DIR=${LIBRARY_DIR}/petsc
    PETSC_ARCH=debug

    git clone -b release https://gitlab.com/petsc/petsc.git ${PETSC_DIR}
    cd ${PETSC_DIR}

    ./configure PETSC_ARCH=debug --with-scalapack-include=$MKLROOT/include --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" --with-mpi-dir=$I_MPI_ROOT --with-mkl=$MKLROOT --with-blas-lapack-dir=$MKLROOT --with-mkl_pardiso-dir=$MKLROOT --with-mkl_cpardiso-dir=$MKLROOT --download-hdf5 --download-mumps --with-clanguage=cxx --with-debugging=1
    make all check

    PETSC_ARCH=optimised
    ./configure PETSC_ARCH=optimised --with-scalapack-include=$MKLROOT/include --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" --with-mpi-dir=$I_MPI_ROOT --with-mkl=$MKLROOT --with-blas-lapack-dir=$MKLROOT --with-mkl_pardiso-dir=$MKLROOT --with-mkl_cpardiso-dir=$MKLROOT --download-hdf5 --download-mumps --with-clanguage=cxx --with-debugging=0 --COPTFLAGS='-O3' --FOPTFLAGS='-O3' --CXXOPTFLAGS='-O3'
    make all check

    # add to bashrc
    if ! grep -q "^export PETSC_DIR=${PETSC_DIR}" "$rcFile"; then
        echo "export PETSC_DIR=${PETSC_DIR}" >> "$rcFile"
        echo "[inserted] export PETSC_DIR=${PETSC_DIR}"
    fi
    if ! grep -q "^export PETSC_ARCH=${PETSC_ARCH}" "$rcFile"; then
        echo "export PETSC_ARCH=${PETSC_ARCH}" >> "$rcFile"
        echo "[inserted] export PETSC_ARCH=${PETSC_ARCH}"
    fi

    cd ${BASEFOLDER}
fi

#rapidjson
RAPIDJSON_DIR=${LIBRARY_DIR}/rapidjson
if [ ! -d $RAPIDJSON_DIR ]; then
    git clone https://github.com/Tencent/rapidjson.git ${RAPIDJSON_DIR}
fi

#eigen3
EIGEN3_Dir=${LIBRARY_DIR}/eigen
if [ ! -d $EIGEN3_Dir ]; then
    git clone https://gitlab.com/libeigen/eigen.git ${EIGEN3_Dir}
fi

#highfive
HIGHFIVE_DIR=${LIBRARY_DIR}/HighFive
if [ ! -d $HIGHFIVE_DIR ]; then
    git clone https://github.com/BlueBrain/HighFive.git ${HIGHFIVE_DIR}
fi

# #VTK
VTK_DIR=${LIBRARY_DIR}/VTK
if [ ! -d $VTK_DIR ]; then
    git clone --recursive https://gitlab.kitware.com/vtk/vtk.git ${VTK_DIR}
    cd ${VTK_DIR}
    mkdir build
    cd build
    cmake -GNinja ../
    cmake --build .
    cd ${BASEFOLDER}

    if ! grep -q "^export VTK_DIR=${VTK_DIR}" "$rcFile"; then
        echo "export VTK_DIR=${VTK_DIR}" >> "$rcFile"
        echo "[inserted] export VTK_DIR=${VTK_DIR}"
    fi
fi

# if ! grep -q "^export LIBGL_ALWAYS_INDIRECT=0" "$rcFile"; then
#     echo "export LIBGL_ALWAYS_INDIRECT=0" >> "$rcFile"
#     echo "[inserted] LIBGL_ALWAYS_INDIRECT=0"
# fi
