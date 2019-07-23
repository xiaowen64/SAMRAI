set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/thirdparty/")

# MPI is setup by BLT
if (MPI_FOUND)
  set(HAVE_MPI True)
else ()
  set(LACKS_MPI True)
endif ()

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


#HAVE_HYPRE
if (ENABLE_HYPRE)
  find_package(HYPRE REQUIRED)

  if(HYPRE_FOUND)
    set (HAVE_HYPRE True)

    blt_register_library(
      NAME HYPRE
      INCLUDES ${HYPRE_INCLUDE_DIRS}
      LIBRARIES ${HYPRE_LIBRARIES})
  endif ()
endif ()

# OpenMP
if (ENABLE_OPENMP)
  if (OPENMP_FOUND)
    set(HAVE_OPENMP True)
  endif ()
endif ()

#HAVE_PETSC
if (ENABLE_PETSC)
  find_package(PETSc REQUIRED)

  if (PETSC_FOUND)
    set (HAVE_PETSC True)

    blt_register_library(
      NAME PETSc
      INCLUDES ${PETSC_INCLUDE_DIRS}
      LIBRARIES ${PETSC_LIBRARIES})
  endif ()
endif()

#HAVE_SILO
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

#HAVE_SUNDIALS
if (ENABLE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
  if (SUNDIALS_FOUND)
    set (HAVE_SUNDIALS True)

    blt_register_library(
      NAME SUNDIALS
      INCLUDES ${SUNDIALS_INCLUDE_DIRS}
      LIBRARIES ${SUNDIALS_LIBRARIES})
  endif ()
endif ()

#HAVE_CONDUIT
if (ENABLE_CONDUIT)
  find_package(CONDUIT REQUIRED)
  if (CONDUIT_FOUND)
    set (HAVE_CONDUIT True)
  endif ()
endif ()
