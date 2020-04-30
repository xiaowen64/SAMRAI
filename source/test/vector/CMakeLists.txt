if(ENABLE_MPI)
  set(TASKS 1)
else()
  set(TASKS 0)
endif()

blt_add_executable(
  NAME kvtest
  SOURCES kvtest.C
  DEPENDS_ON
    SAMRAI_solv
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_hier
    SAMRAI_tbox)

blt_add_test(
  NAME kvtest
  COMMAND kvtest
  NUM_MPI_TASKS ${TASKS})

blt_add_executable(
  NAME pvtest
  SOURCES pvtest.C
  DEPENDS_ON
    SAMRAI_solv
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_hier
    SAMRAI_tbox)

blt_add_test(
  NAME pvtest
  COMMAND pvtest
  NUM_MPI_TASKS ${TASKS})

