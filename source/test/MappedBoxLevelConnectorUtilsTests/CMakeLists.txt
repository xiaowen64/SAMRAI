set (mapped_boxlevel_sources
  main-mblcu.C)

blt_add_executable(
  NAME mblcu
  SOURCES ${mapped_boxlevel_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

target_compile_definitions(mblcu PUBLIC TESTING=1)
