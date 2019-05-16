set (mblkeuler_sources
  main.C
  MblkEuler.C
  MblkGeometry.C
  boundaryconditions.f)

set (mblkeuler_depends_on
    testlib
    SAMRAI_appu
    SAMRAI_algs
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

blt_add_executable(
  NAME mblkeuler
  SOURCES ${mblkeuler_sources}
  DEPENDS_ON
    ${mblkeuler_depends_on})

target_compile_definitions(mblkeuler PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME mblkeuler
  EXECUTABLE mblkeuler
  INPUTS ${test_inputs}
  PARALLEL TRUE)

