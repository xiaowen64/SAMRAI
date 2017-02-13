set (indexdata_sources
  main.C)

blt_add_executable(
  NAME indexdata
  SOURCES ${indexdata_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_pdat)

target_compile_definitions(indexdata PUBLIC TESTING=1)
