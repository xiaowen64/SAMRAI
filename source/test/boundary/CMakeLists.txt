set ( boundary_sources
  main.C
  BoundaryDataTester.C)

blt_add_executable(
  NAME boundary
  SOURCES ${boundary_sources}
  DEPENDS_ON
    SAMRAI_tbox)

target_include_directories( boundary
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/boundary)

target_compile_definitions(boundary TESTING=1)
