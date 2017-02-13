set ( fill_pattern_sources
  main-fill-pattern.C)

blt_add_executable(
  NAME fill_pattern
  SOURCES ${fill_pattern_sources}
  DEPENDS_ON
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_xfer
    SAMRAI_hier
    SAMRAI_tbox)

target_compile_definitions(fill_pattern PUBLIC TESTING=1)
