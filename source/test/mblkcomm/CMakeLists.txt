set (mblkcomm_sources
  CellMultiblockTest.C
  EdgeMultiblockTest.C
  FaceMultiblockTest.C
  main.C
  MultiblockTester.C
  NodeMultiblockTest.C
  PatchMultiblockTestStrategy.C
  SideMultiblockTest.C)

set (mblkcomm_depends_on
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_mesh
  SAMRAI_tbox)

blt_add_executable(
  NAME mblkcomm
  SOURCES ${mblkcomm_sources}
  DEPENDS_ON
    ${mblkcomm_depends_on})

target_compile_definitions(mblkcomm PUBLIC TESTING=1)
