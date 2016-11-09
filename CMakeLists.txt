cmake_minimum_required(VERSION 3.1)

project(SAMRAI C CXX Fortran)

set(SAMRAI_VERSION_MAJOR 3)
set(SAMRAI_VERSION_MINOR 11)
set(SAMRAI_VERSION_PATCHLEVEL 1)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "The BLT submodule is not present. \
If in the git repository, please run the following commands:\n \
git submodule init\n \
git submodule update")
endif ()

# Fortran enabled for BLT setup
option(ENABLE_FORTRAN "Enable Fortran" On)

option(ENABLE_OPENMP "Enable OpenMP" On)

# Third-party packages
option(ENABLE_MPI "Enable MPI." On)
option(ENABLE_HDF5 "Enable HDF5." On)
option(ENABLE_HYPRE "Enable HYPRE" Off)
option(ENABLE_PETSC "Enable PETSc" Off)
option(ENABLE_PTSCOTCH "Enable Scotch" Off)
option(ENABLE_SILO "Enable Silo" Off)
option(ENABLE_SUNDIALS "Enable SUNDIALS" Off)

# SAMRAI options
option(ENABLE_CHECK_ASSERTIONS "Enable assertion checking." Off)
option(ENABLE_CHECK_DEV_ASSERTIONS "Enable SAMRAI developer assertion checking." Off)
option(ENABLE_CHECK_DIM_ASSERTIONS "Enable assertion checking for dimensions." Off)
option(ENABLE_BOX_COUNTING "Turns on box telemetry." Off)
option(ENABLE_DEPRECATED "Build with deprecated features." On)
option(ENABLE_TIMERS "Enable SAMRAI timers." On)
option(DEBUG_INITIALIZE_UNDEFINED "Initialize new memory to undefined values" Off)
set(MAXDIM 3 CACHE INT "Maximum allowed dimension")

include(blt/SetupBLT.cmake)

include(cmake/CMakeBasics.cmake)

add_subdirectory(source)

add_subdirectory(tools)

add_subdirectory(docs)
