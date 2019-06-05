set (hierarchy_sources
  main.C
  HierarchyTester.C)

blt_add_executable(
  NAME hierarchy
  SOURCES ${hierarchy_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_tbox
    SAMRAI_mesh
    SAMRAI_geom)

target_compile_definitions(hierarchy PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME hierarchy
  EXECUTABLE hierarchy
  INPUTS ${test_inputs}
  PARALLEL TRUE)
