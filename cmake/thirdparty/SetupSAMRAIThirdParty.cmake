set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/thirdparty/")

if (MPI_FOUND)
  set(HAVE_MPI True)
else ()
  set(LACKS_MPI True)
endif ()

find_package(Boost REQUIRED)
if (${Boost_FOUND})
  set(HAVE_BOOST_HEADERS True)
endif ()

blt_register_library(
  NAME boost
  INCLUDES ${Boost_INCLUDE_DIR})

if (ENABLE_HDF5)
  if (NOT ENABLE_MPI)
    message(FATAL_ERROR "HDF5 requires MPI.")
  endif ()

  find_package(HDF5)

  if(HDF5_FOUND)
    set (HAVE_HDF5 True)

    blt_register_library(
      NAME hdf5
      INCLUDES ${HDF5_INCLUDE_DIRS}
      LIBRARIES ${HDF5_LIBRARIES})
  endif ()
endif ()

# find_package(

#HAVE_HYPRE
if (ENABLE_HYPRE)
  find_package(HYPRE)

  if(HYPRE_FOUND)
    set (HAVE_HPYRE True)

    blt_register_library(
      NAME hypre
      INCLUDES ${HYPRE_INCLUDE_DIRS}
      LIBRARIES ${HPYRE_LIBRARIES})
  endif ()
endif ()

#HAVE_OPENMP
if (ENABLE_OPENMP)
  if (OPENMP_FOUND)
    set(HAVE_OPENMP True)
  endif ()
endif ()

#HAVE_PETSC
if (ENABLE_PETSC)
  find_package(PETSc)

  if (PETSC_FOUND)
    set (HAVE_PETSC True)

    blt_register_library(
      NAME PETSc
      INCLUDES ${PETSC_INCLUDES}
      LIBRARIES ${PETSC_LIBRARIES})
  endif ()
endif()

#HAVE_PTSCOTCH
if (ENABLE_PTSCOTCH)
endif ()

#HAVE_SILO
if (ENABLE_SILO)
  find_package(silo)

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
endif ()


#HAVE_TAU
if (ENABLE_TAU)
endif ()

#HAVE_VAMPIR
if (ENABLE_VAMPIR)
endif ()
