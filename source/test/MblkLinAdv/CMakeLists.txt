set (mblklinadv_sources
  MblkGeometry.C
  MblkLinAdv.C
  SkeletonOutersideDoubleWeightedAverage.C
  main.C)

set (mblklinadv_depends_on
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_pdat
    SAMRAI_algs
    SAMRAI_appu
    SAMRAI_tbox)

blt_add_executable(
  NAME mblklinadv
  SOURCES ${mblklinadv_sources}
  DEPENDS_ON
    ${mblklinadv_depends_on})

target_compile_definitions(mblklinadv PUBLIC TESTING=1)
