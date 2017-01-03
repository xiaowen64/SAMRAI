set (hierarchy_sources
  main.C
  HierarchyTester.C)

blt_add_executable(
  NAME hierarchy
  SOURCES ${hierarchy_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_tbox
    SAMRAI_mesh
    SAMRAI_geom)

target_compile_definitions(hierarchy PUBLIC TESTING=1)
