set ( async_br_sources
  main.C
  ABRTest.C)

blt_add_executable(
  NAME async_br
  SOURCES ${async_br_sources}
  DEPENDS_ON
    testlib
    SAMRAI_appu
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

target_include_directories( async_br
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/clustering/async_br)

target_compile_definitions(async_br PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME async_br
  EXECUTABLE async_br
  INPUTS ${test_inputs}
  PARALLEL TRUE)
