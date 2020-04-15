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

set(communication_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(communication_depends ${communication_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(communication_depends ${communication_depends} cuda)
endif ()

set(communication_depends ${communication_depends} testlib)


blt_add_executable(
  NAME communication
  SOURCES ${communication_sources}
  DEPENDS_ON ${communication_depends})

target_include_directories( communication
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/communication)

target_compile_definitions(communication PUBLIC TESTING=1)

set_target_properties(communication PROPERTIES LINKER_LANGUAGE CXX)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME communication
  EXECUTABLE communication
  INPUTS ${test_inputs}
  PARALLEL TRUE)
