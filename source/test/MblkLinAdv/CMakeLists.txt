set (mblklinadv_sources
  MblkGeometry.C
  MblkLinAdv.C
  SkeletonOutersideDoubleWeightedAverage.C
  main.C)

set (mblklinadv_depends_on
    testlib
    SAMRAI_appu
    SAMRAI_algs
    SAMRAI_solv
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_math
    SAMRAI_pdat
    SAMRAI_xfer
    SAMRAI_hier
    SAMRAI_tbox)

blt_add_executable(
  NAME mblklinadv
  SOURCES ${mblklinadv_sources}
  DEPENDS_ON
    ${mblklinadv_depends_on})

  #target_compile_definitions(mblklinadv PUBLIC TESTING=1)
