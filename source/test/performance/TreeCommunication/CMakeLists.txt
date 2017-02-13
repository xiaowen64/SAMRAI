set (treecommunication_sources
  tcb.C)

blt_add_executable(
  NAME treecommunication
  SOURCES ${treecommunication_sources}
  DEPENDS_ON
    SAMRAI_tbox)

target_compile_definitions(treecommunication PUBLIC TESTING=1)
