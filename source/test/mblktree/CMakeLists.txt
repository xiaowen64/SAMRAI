set (mblktree_sources
  main-mbtree.C)

set (mblktree_depends_on
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox)

blt_add_executable(
  NAME mblktree
  SOURCES ${mblktree_sources}
  DEPENDS_ON
    ${mblktree_depends_on})

target_compile_definitions(mblktree PUBLIC TESTING=1)
