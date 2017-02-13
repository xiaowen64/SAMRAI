set ( communication_sources
  main.C
  CommTester.C
  PatchDataTestStrategy.C
  CellDataTest.C
  EdgeDataTest.C
  FaceDataTest.C
  NodeDataTest.C
  OuterfaceDataTest.C
  OuternodeDataTest.C
  OutersideDataTest.C
  SideDataTest.C)

blt_add_executable(
  NAME communication
  SOURCES ${communication_sources}
  DEPENDS_ON
    testlib
    SAMRAI_appu
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_mesh
    SAMRAI_xfer
    SAMRAI_geom)

target_include_directories( communication
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/communication)

target_compile_definitions(communication PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME communication
  EXECUTABLE communication
  INPUTS ${test_inputs}
  PARALLEL True)
