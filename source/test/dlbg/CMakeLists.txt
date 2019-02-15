set ( dlbg_sources
  main.C
  DLBGTest.C)

blt_add_executable(
  NAME dlbg
  SOURCES ${dlbg_sources}
  DEPENDS_ON
    testlib
    SAMRAI_tbox
    SAMRAI_appu
    SAMRAI_geom
    SAMRAI_hier
    SAMRAI_solv
    SAMRAI_math
    SAMRAI_mesh)

target_compile_definitions(dlbg PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME dlbg
  EXECUTABLE dlbg
  INPUTS ${test_inputs}
  PARALLEL TRUE)
