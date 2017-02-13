set (meshgeneration_sources
  mg.C)

blt_add_executable(
  NAME meshgeneration
  SOURCES ${meshgeneration_sources}
  DEPENDS_ON
    SAMRAI)

target_compile_definitions(meshgeneration PUBLIC TESTING=1)
