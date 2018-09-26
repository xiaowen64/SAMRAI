cmake_minimum_required(VERSION 3.1)

project(SAMRAI C CXX Fortran)

set(SAMRAI_VERSION_MAJOR 3)
set(SAMRAI_VERSION_MINOR 11)
set(SAMRAI_VERSION_PATCHLEVEL 1)

cmake_policy(SET CMP0057 NEW)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "The BLT submodule is not present. \
If in the git repository, please run the following commands:\n \
git submodule init\n \
git submodule update")
endif ()

# BLT Options need to be set in the CACHE
set(ENABLE_FORTRAN On CACHE Bool "Enable Fortran")
set(ENABLE_OPENMP On CACHE Bool "Enable OpenMP")
set(ENABLE_CUDA On CACHE Bool "Enable CUDA")
set(ENABLE_MPI On CACHE Bool "Enable MPI.")
set(ENABLE_COPY_HEADERS On CACHE Bool "Copy headers")

# Third party libraries
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

# Converts C-style comments to Fortran style in preprocessed files
set(BLT_FORTRAN_FLAGS "-WF,-C!" CACHE STRING "")

include_directories(${CMAKE_BINARY_DIR}/include)

include(cmake/CMakeBasics.cmake)

add_subdirectory(source)
add_subdirectory(tools)
add_subdirectory(docs)

install(DIRECTORY
  ${CMAKE_BINARY_DIR}/include
  DESTINATION ${CMAKE_INSTALL_PREFIX})

add_custom_target(SAMRAI)
add_dependencies(SAMRAI
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_geom
  SAMRAI_hier
  SAMRAI_math
  SAMRAI_mesh
  SAMRAI_pdat
  SAMRAI_solv
  SAMRAI_tbox
  SAMRAI_xfer)
