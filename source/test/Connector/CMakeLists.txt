set ( connector_sources
  main.C)

blt_add_executable(
  NAME Connector
  SOURCES ${connector_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_geom)

target_include_directories( Connector
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/Connector)

target_compile_definitions(Connector PUBLIC TESTING=1)
