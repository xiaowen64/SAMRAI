set(ENABLE_MPI On CACHE BOOL "")
set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-7.1.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-7.1.0/bin/mpic++" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-7.1.0/bin/mpif90" CACHE PATH "")

set(ENABLE_HDF5 Off CACHE BOOL "")
