set (mapping_connector_sources
  main-modify.C)

blt_add_executable(
  NAME mapping-connector
  SOURCES ${mapping_connector_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

target_compile_definitions(mapping-connector PUBLIC TESTING=1)
