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
  COMMAND kvtest)

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
  COMMAND pvtest)
