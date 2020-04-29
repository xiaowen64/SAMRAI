set (indexdata_sources
  main.C)

blt_add_executable(
  NAME indexdata
  SOURCES ${indexdata_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_pdat)

target_compile_definitions(indexdata PUBLIC TESTING=1)

if(ENABLE_MPI)
  set(TASKS 1)
else()
  set(TASKS 0)
endif()

blt_add_test(
  NAME indexdata
  COMMAND indexdata
  NUM_MPI_TASKS ${TASKS})
