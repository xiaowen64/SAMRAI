set(CMAKE_Fortran_FORMAT FIXED)

# set (CUDA_NVCC_FLAGS
#   -Xcompiler; ${OpenMP_CXX_FLAGS};
#   -g; -G; -arch sm_60
#   --expt-extended-lambda;
#   -std=c++11)

set (CUDA_NVCC_FLAGS
  -Xcompiler; ${OpenMP_CXX_FLAGS};
  -O3;
  -arch sm_60;
  --expt-extended-lambda;
  -std=c++11)
