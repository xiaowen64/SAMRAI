set (mblkeuler_sources
  main.C
  MblkEuler.C
  MblkGeometry.C
  bc_common.f
  boundaryconditions.f)

set (mblkeuler_depends_on
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

blt_add_executable(
  NAME mblkeuler
  SOURCES ${mblkeuler_sources}
  DEPENDS_ON
    ${mblkeuler_depends_on})

target_compile_definitions(mblkeuler PUBLIC TESTING=1)
