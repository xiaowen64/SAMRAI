blt_add_executable(
  NAME transformation
  SOURCES main.C
  DEPENDS_ON
    SAMRAI_pdat
    SAMRAI_hier
    SAMRAI_tbox)

target_include_directories( transformation
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/transformation)
