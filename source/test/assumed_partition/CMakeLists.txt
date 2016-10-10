blt_add_executable(
  NAME assumed_partition
  SOURCES main.C
  DEPENDS_ON
    SAMRAI_geom
    SAMRAI_hier
    SAMRAI_tbox)

target_include_directories( assumed_partition
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/assumed_partition)
