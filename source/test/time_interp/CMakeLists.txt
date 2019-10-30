set ( time_interp_sources
  main.C
)

blt_add_executable(
  NAME time_interp
  SOURCES ${time_interp_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_pdat)

target_include_directories( time_interp
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/time_interp)

target_compile_definitions(time_interp PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME time_interp
  EXECUTABLE time_interp
  INPUTS ${test_inputs}
  PARALLEL FALSE)
