set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/thirdparty/")

if (NOT ENABLE_UMPIRE)
  message(FATAL_ERROR "This code branch requires UMPIRE, for now...")
endif ()

if (ENABLE_CUDA)
  if (NOT ENABLE_RAJA)
    message(FATAL_ERROR "CUDA support requires RAJA")
  endif ()

  if (NOT ENABLE_UMPIRE)
    message(FATAL_ERROR "CUDA support requires UMPIRE")
  endif ()
endif ()

# MPI is setup by BLT
if (MPI_FOUND)
  set(HAVE_MPI True)
else ()
  set(LACKS_MPI True)
endif ()

# OpenMP is set up by BLT
if (ENABLE_OPENMP)
  if (OPENMP_FOUND)
    set(HAVE_OPENMP True)
  endif ()
endif ()

# CUDA is setup by BLT
if (ENABLE_CUDA)
  if (CUDA_FOUND)
    set (HAVE_CUDA True)
  endif ()
endif ()

# UMPIRE
if (ENABLE_UMPIRE)
  find_package(umpire REQUIRED)

  set (HAVE_UMPIRE True)

  blt_register_library(
    NAME umpire
    INCLUDES ${UMPIRE_INCLUDE_DIRS}
    LIBRARIES umpire)
endif ()

# RAJA
if (ENABLE_RAJA)
  if (NOT ENABLE_UMPIRE)
    message(FATAL_ERROR "RAJA support requires UMPIRE.")
  endif ()

  find_package(RAJA REQUIRED)

  set (raja_depends_on)
  if (ENABLE_CUDA)
    list (APPEND raja_depends cuda)
  endif ()

  if (ENABLE_OPENMP)
    list (APPEND raja_depends openmp)
  endif ()

  if (RAJA_FOUND)
    set (HAVE_RAJA True)

    blt_register_library(
      NAME RAJA
      INCLUDES ${RAJA_INCLUDE_DIR}
      LIBRARIES RAJA
      DEPENDS_ON ${raja_depends})
  endif ()
endif ()

# HDF5
if (ENABLE_HDF5)
  if (NOT ENABLE_MPI)
    message(FATAL_ERROR "HDF5 requires MPI.")
  endif ()

  find_package(HDF5 REQUIRED)

  if(HDF5_FOUND)
    set (HAVE_HDF5 True)

    blt_register_library(
      NAME hdf5
      INCLUDES ${HDF5_INCLUDE_DIRS}
      LIBRARIES ${HDF5_C_LIBRARIES})
  endif ()
endif ()

# HYPRE
if (ENABLE_HYPRE)
  find_package(HYPRE REQUIRED)
  # TODO: Ensure this is set in SAMRAI_config.h...

  if(HYPRE_FOUND)
    set (HAVE_HPYRE True)

    blt_register_library(
      NAME hypre
      INCLUDES ${HYPRE_INCLUDE_DIRS}
      LIBRARIES ${HYPRE_LIBRARIES})
  endif ()
endif ()

# PETSC
if (ENABLE_PETSC)
  find_package(PETSc REQUIRED)

  if (PETSC_FOUND)
    set (HAVE_PETSC True)

    blt_register_library(
      NAME PETSc
      INCLUDES ${PETSC_INCLUDES}
      LIBRARIES ${PETSC_LIBRARIES})
  endif ()
endif()

# PTSCOTCH
if (ENABLE_PTSCOTCH)
  find_package(Scotch REQUIRED)

  if (Scotch_FOUND)
    set (HAVE_SCOTCH True)

    blt_register_library(
      NAME Scotch
      INCLUDES ${SCOTCH_INCLUDES}
      LIBRARIES ${SCOTCH_LIBRARIES})
  endif ()
endif ()

# SILO
if (ENABLE_SILO)
  find_package(SILO REQUIRED)

  if (SILO_FOUND)
    set (HAVE_SILO True)

    blt_register_library(
      NAME silo
      INCLUDES ${SILO_INCLUDE_DIRS}
      LIBRARIES ${SILO_LIBRARIES})
  endif ()
endif ()

# SUNDIALS
if (ENABLE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
  if (SUNDIALS_FOUND)
    set (HAVE_SUNDIALS True)

    blt_register_library(
      NAME SUNDIALS
      INCLUDES ${SUNDIALS_INCLUDES}
      LIBRARIES ${SUNDIALS_LIBRARIES})
  endif ()
endif ()

# CONDUIT
if (ENABLE_CONDUIT)
  find_package(CONDUIT REQUIRED)
  if (CONDUIT_FOUND)
    set (HAVE_CONDUIT True)
  endif ()
endif ()
