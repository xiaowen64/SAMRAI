if (${MPI_FOUND})
  set(HAVE_MPI True)
else ()
  set(LACKS_MPI True)
endif ()

find_package(Boost REQUIRED)

blt_register_library(
  NAME boost
  INCLUDES ${Boost_INCLUDE_DIR})

if (${ENABLE_HDF5})
  if (NOT ${ENABLE_MPI})
    message(FATAL_ERROR "HDF5 requires MPI.")
  endif ()

  find_package(HDF5)

  if(${HDF5_FOUND})
    set (HAVE_HDF5 True)

    blt_register_library(
      NAME hdf5
      INCLUDES ${HDF5_INCLUDE_DIRS}
      LIBRARIES ${HDF5_LIBRARIES})
  endif ()
endif ()
