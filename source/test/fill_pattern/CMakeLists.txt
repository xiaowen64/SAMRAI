set ( fill_pattern_sources
  main-fill-pattern.C)

set(fill_pattern_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(fill_pattern_depends ${fill_pattern_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(fill_pattern_depends ${fill_pattern_depends} cuda)
endif ()


blt_add_executable(
  NAME fill_pattern
  SOURCES ${fill_pattern_sources}
  DEPENDS_ON ${fill_pattern_depends})

target_compile_definitions(fill_pattern PUBLIC TESTING=1)

blt_add_test(
  NAME fill_pattern
  COMMAND fill_pattern)

