set (indexdata_sources
  main.C)

set(indexdata_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(indexdata_depends ${indexdata_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(indexdata_depends ${indexdata_depends} cuda)
endif ()

blt_add_executable(
  NAME indexdata
  SOURCES ${indexdata_sources}
  DEPENDS_ON ${indexdata_depends})

target_compile_definitions(indexdata PUBLIC TESTING=1)

blt_add_test(
  NAME indexdata
  COMMAND indexdata)
