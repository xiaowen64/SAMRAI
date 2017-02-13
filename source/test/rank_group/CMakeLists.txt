set (rank_group_sources
  main.C)

blt_add_executable(
  NAME rank_group
  SOURCES ${rank_group_sources}
  DEPENDS_ON
  SAMRAI_mesh
  SAMRAI_hier
  SAMRAI_geom
  SAMRAI_pdat
  SAMRAI_appu
  SAMRAI_tbox)

target_compile_definitions(rank_group PUBLIC TESTING=1)
