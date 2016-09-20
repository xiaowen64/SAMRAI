blt_add_executable(NAME restart-redistribute
  SOURCES main.C RedistributedRestartUtility.C
  DEPENDS_ON SAMRAI_tbox)

target_include_directories(
  restart-redistribute
  PRIVATE ${PROJECT_SOURCE_DIR}/source/)
