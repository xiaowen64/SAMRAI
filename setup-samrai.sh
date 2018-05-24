#!/usr/bin/env bash

arch=$(uname -m)
compiler=clang
scratch_dir=${HOME}/Projects/SAMRAI/
build_prefix=${scratch_dir}/build/${arch}-${compiler}/
install_prefix=${scratch_dir}/install/${arch}-${compiler}/
source_prefix=${HOME}/Projects/SAMRAI/code

module load spectrum-mpi/2018.02.05
module load cuda/9.0.184

pkg=$1
if [[ "$1" == "-c" ]]; then
    pkg=$2
    dir=${build_prefix}/$pkg
    echo Removing $dir ...
    rm -rf $dir
fi

if [[ "$compiler" == "xl" ]]; then
    cc=xlc
    cpp=xlC
    fc=xlf
    extra_cflags=""
elif [[ "$compiler" == "clang" ]]; then
    cc=clang
    cpp=clang++
    fc=xlf
    extra_cflags=""
else
    exit 1
fi

build_dir=${build_prefix}/${pkg}
install_dir=${install_prefix}/${pkg}
if [ "$pkg" == "RAJA" ]; then
    source_dir=${source_prefix}/${pkg}
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DENABLE_OPENMP=Off \
          -DENABLE_CUDA=ON \
          -DENABLE_TESTS=OFF \
          -DENABLE_NESTED=ON \
          -DRAJA_CUDA_ARCH=sm_60 \
          $source_dir
elif [ "$pkg" == "umpire" ]; then
    source_dir=${source_prefix}/${pkg}
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_BUILD_TYPE=Release \
          -DENABLE_OPENMP=On \
          -DENABLE_CUDA=ON \
          -DOpenMP_Fortran_FLAGS="-qsmp" \
          -DOpenMP_Fortran_LIB_NAMES="" \
          -DBUILD_SHARED_LIBS=OFF \
          -DENABLE_FORTRAN=Off \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          $source_dir
elif [ "$pkg" == "HDF5" ]; then
    source_dir=${source_prefix}/hdf5-1.10.1
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Release \
          -DHDF5_BUILD_EXAMPLES=OFF \
          $source_dir
elif [ "$pkg" == "samrai-uvm" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir && {
        branch=$(git rev-parse --abbrev-ref HEAD)
        [[ "$branch" != "feature/uvm" ]] && exit 1
    }
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_Fortran_COMPILER=$fc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_PREFIX_PATH="${install_prefix}/RAJA/share/raja/cmake;${install_prefix}/HDF5/share/cmake" \
          -Dumpire_DIR="${install_prefix}/umpire/share/umpire/cmake" \
          -DMPI_Fortran_COMPILER=mpixlf \
          -DENABLE_HDF5=On \
          -DENABLE_MPI=On \
          -DENABLE_OPENMP=ON \
          -DOpenMP_Fortran_FLAGS="-qsmp" \
          -DOpenMP_Fortran_LIB_NAMES="" \
          -DENABLE_CUDA=ON \
          -DBUILD_SHARED_LIBS=OFF \
          -DENABLE_GTEST=ON \
          -DBUILD_GTEST=ON \
          -DBUILD_TESTING=ON \
          -DENABLE_EXAMPLES=OFF \
          $source_dir
fi
