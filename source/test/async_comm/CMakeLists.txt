blt_add_executable(
  NAME main-async_comm
  SOURCES main-async_comm.C
  DEPENDS_ON
    SAMRAI_tbox)

blt_add_executable(
  NAME main-peer_comm
  SOURCES main-peer_comm.C
  DEPENDS_ON
    SAMRAI_tbox)

target_include_directories( main-async_comm
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/main-async_comm)

target_include_directories( main-peer_comm
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/main-peer_comm)
