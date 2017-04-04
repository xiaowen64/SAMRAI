set (meshgeneration_sources
  mg.C)

set (meshgeneration_depends_on
  testlib
  SAMRAI_mesh
  SAMRAI_hier
  SAMRAI_geom
  SAMRAI_pdat
  SAMRAI_appu
  SAMRAI_tbox)

blt_add_executable(
  NAME meshgeneration
  SOURCES ${meshgeneration_sources}
  DEPENDS_ON ${meshgeneration_depends_on})

target_compile_definitions(meshgeneration PUBLIC TESTING=1)
