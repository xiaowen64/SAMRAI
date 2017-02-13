set ( boundary_sources
  main.C
  BoundaryDataTester.C)

blt_add_executable(
  NAME boundary
  SOURCES ${boundary_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_appu
    SAMRAI_xfer
    SAMRAI_geom)

target_include_directories( boundary
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/boundary)

target_compile_definitions(boundary PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME boundary
  EXECUTABLE boundary
  INPUTS ${test_inputs}
  PARALLEL True)
