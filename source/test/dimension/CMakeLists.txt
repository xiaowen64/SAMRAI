set ( dimension_sources
  main.C)

blt_add_executable(
  NAME dimension
  SOURCES ${dimension_sources}
  DEPENDS_ON
    SAMRAI_tbox)

target_compile_definitions(dimension PUBLIC TESTING=1)

blt_add_test(
  NAME dimension
  COMMAND dimension)
