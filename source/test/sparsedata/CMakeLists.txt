set (sparsedata_sources
  main.C
  SparseDataTester.C)

blt_add_executable(
  NAME main-sparsedata
  SOURCES ${sparsedata_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_pdat
    SAMRAI_geom)

target_compile_definitions(main-sparsedata PUBLIC TESTING=1)

set (sparsedataops_sources
  sparse_dataops.C)

blt_add_executable(
  NAME sparse
  SOURCES ${sparsedataops_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_hier
    SAMRAI_pdat
    SAMRAI_geom)

target_compile_definitions(sparse PUBLIC TESTING=1)
