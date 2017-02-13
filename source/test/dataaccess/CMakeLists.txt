set ( dataaccess_sources
  main.C)

blt_add_executable(
  NAME dataaccess
  SOURCES ${dataaccess_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_pdat
    SAMRAI_hier)

target_include_directories( dataaccess
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataaccess)

target_compile_definitions(dataaccess PUBLIC TESTING=1)
