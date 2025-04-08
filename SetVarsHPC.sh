#!/bin/bash

echo Intel
module load intel/2023a

echo package config
module load pkg-config/0.29.2

echo CMake
module load CMake/3.26.3-GCCcore-12.3.0

#echo HDF5
#module load HDF5/1.14.0-iimpi-2023a

#echo RapidJSON
#module load RapidJSON/1.1.0-20230928-GCCcore-12.3.0

#echo HighFive
#module load HighFive

#echo Eigen
#module load Eigen/3.4.0-GCCcore-13.2.0

echo PETSC
export PETSC_DIR=$DATA/Codes/petsc
export PETSC_ARCH=optimised

export CurrentEnv="hpc"
export LIBRARY_DIR=$DATA/Codes/Libraries

echo All Finished Setting up
