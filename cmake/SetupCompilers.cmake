set(CMAKE_Fortran_FORMAT FIXED)

set (CUDA_NVCC_FLAGS
  -Xcompiler; ${OpenMP_CXX_FLAGS};
  -G; -g; -arch sm_60
  --expt-extended-lambda;
  -std=c++11)
