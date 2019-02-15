blt_add_executable(NAME variables
  SOURCES main-var_db.C
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_pdat
    SAMRAI_hier
  )

target_compile_definitions(variables PUBLIC TESTING=1)

blt_add_test(
  NAME variables
  COMMAND variables)
