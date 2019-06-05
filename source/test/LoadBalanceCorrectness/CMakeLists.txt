set(lbcorrectness_sources
  main-lbcorrectness.C)

blt_add_executable(
  NAME lbcorrectness
  SOURCES ${lbcorrectness_sources}
  DEPENDS_ON
    testlib
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_appu
    SAMRAI_mesh
    SAMRAI_pdat
    SAMRAI_geom)

target_compile_definitions(lbcorrectness PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)
file (GLOB test_baselines ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*baseline*)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/test_inputs)
execute_process(COMMAND ln -s ${test_baselines} ${CMAKE_CURRENT_BINARY_DIR}/test_inputs)

samrai_add_tests(
  NAME lbcorrectness
  EXECUTABLE lbcorrectness
  INPUTS ${test_inputs}
  PARALLEL TRUE)
