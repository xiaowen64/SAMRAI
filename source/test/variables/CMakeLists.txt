blt_add_executable(NAME variables
  SOURCES main-var_db.C
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_pdat
    SAMRAI_hier
  )
