#!/usr/bin/env bash

arch=$(uname -m)
compiler=clang
build_type=Debug
scratch_dir=${HOME}/Projects/Cleverleaf/build
# if [ -z ${SCRATCHDIR+x} ]; then
#     scratch_dir=${HOME}/Projects/SAMRAI
# else
#     scratch_dir=${SCRATCHDIR}
# fi
build_prefix=${scratch_dir}/build/${arch}-${compiler}/${build_type}
install_prefix=${scratch_dir}/install/${arch}-${compiler}/${build_type}
source_prefix=${HOME}/Projects/Cleverleaf

# module load spectrum-mpi/2018.02.05
# module load cuda/9.0.184

echo "Build prefix: $build_prefix"
echo "Install prefix: $install_prefix"

pkg=$1
if [[ "$1" == "-c" ]]; then
    pkg=$2
    dir=${build_prefix}/$pkg
    echo Removing $dir ...
    rm -rf $dir
elif [[ "$1" == "-b" ]]; then
    pkg=$2
    echo ${build_prefix}/${pkg}
    exit 0
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
          -DCMAKE_BUILD_TYPE=$build_type \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DENABLE_OPENMP=OFF \
          -DENABLE_CUDA=ON \
          -DENABLE_TESTS=OFF \
          $source_dir
          # -DRAJA_CUDA_ARCH=sm_60 \
elif [ "$pkg" == "Umpire" ]; then
    source_dir=${source_prefix}/${pkg}
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_BUILD_TYPE=$build_type \
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
          -DCMAKE_BUILD_TYPE=$build_type \
          -DHDF5_BUILD_EXAMPLES=OFF \
          $source_dir
elif [ "$pkg" == "SAMRAI" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir && {
        branch=$(git rev-parse --abbrev-ref HEAD)
	echo "Configuring SAMRAI on branch $branch"
    }
    mkdir -p $build_dir && cd $build_dir
    cmake -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_Fortran_COMPILER=$fc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=$build_type \
          -DRAJA_DIR="${install_prefix}/RAJA/share/raja/cmake" \
          -DHDF5_DIR="${install_prefix}/HDF5/share/cmake" \
          -Dumpire_DIR="${install_prefix}/Umpire/share/umpire/cmake" \
          -DCUDA_NVCC_FLAGS_DEBUG="--expt-relaxed-constexpr" \
          -DCUDA_NVCC_FLAGS_RELEASE="--expt-relaxed-constexpr" \
          -DMPI_Fortran_COMPILER=mpixlf \
          -DENABLE_HDF5=ON \
          -DENABLE_MPI=ON \
          -DENABLE_OPENMP=OFF \
          -DOpenMP_Fortran_FLAGS="-qsmp" \
          -DOpenMP_Fortran_LIB_NAMES="" \
          -DENABLE_CUDA=ON \
          -DBUILD_SHARED_LIBS=OFF \
          -DENABLE_GTEST=ON \
          -DBUILD_GTEST=ON \
          -DBUILD_TESTING=ON \
          -DENABLE_EXAMPLES=OFF \
          $source_dir
elif [ "$pkg" == "cleverleaf" ]; then
    source_dir=${source_prefix}/${pkg}
    cd $source_dir && {
        branch=$(git rev-parse --abbrev-ref HEAD)
	echo "Configuring cleverleaf on branch $branch"
    }
    mkdir -p $build_dir && cd $build_dir
    echo $install_prefix
    cmake \
          -DCMAKE_CXX_COMPILER=$cpp \
          -DCMAKE_C_COMPILER=$cc \
          -DCMAKE_INSTALL_PREFIX=$install_dir \
          -DCMAKE_BUILD_TYPE=$build_type \
          -DCMAKE_CUDA_FLAGS="-arch=sm_70" \
          -DRAJA_DIR="${install_prefix}/RAJA/share/raja/cmake" \
          -DHDF5_DIR="${install_prefix}/HDF5/share/cmake" \
          -Dumpire_DIR="${install_prefix}/Umpire/share/umpire/cmake" \
          -DCMAKE_PREFIX_PATH="${install_prefix}/SAMRAI" \
          -DENABLE_OPENMP=OFF \
          $source_dir
fi
