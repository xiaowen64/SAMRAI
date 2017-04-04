set (oca_sources
  main.C)

set (oca_depends_on
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_tbox)

blt_add_executable(
  NAME overlap-connector-algorithm
  SOURCES ${oca_sources}
  DEPENDS_ON
    ${oca_depends_on})

target_compile_definitions(overlap-connector-algorithm PUBLIC TESTING=1)
