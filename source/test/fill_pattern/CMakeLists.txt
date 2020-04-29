set ( fill_pattern_sources
  main-fill-pattern.C)

blt_add_executable(
  NAME fill_pattern
  SOURCES ${fill_pattern_sources}
  DEPENDS_ON
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_xfer
    SAMRAI_hier
    SAMRAI_tbox)

target_compile_definitions(fill_pattern PUBLIC TESTING=1)

if(ENABLE_MPI)
  set(TASKS 1)
else()
  set(TASKS 0)
endif()

blt_add_test(
  NAME fill_pattern
  COMMAND fill_pattern
  NUM_MPI_TASKS ${TASKS})

