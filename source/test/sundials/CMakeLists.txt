set (sundials_sources
  main.C
  CVODEModel.C)

blt_add_executable(
  NAME sundials_test
  SOURCES ${sundials_sources}
  DEPENDS_ON
    SAMRAI_tbox
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_hier
    SAMRAI_solv)

target_include_directories( sundials_test
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/sundials)

target_compile_definitions(sundials_test PUBLIC TESTING=1)
