set ( async_br_sources
  main.C
  ABRTest.C)

blt_add_executable(
  NAME async_br
  SOURCES ${async_br_sources}
  DEPENDS_ON
  SAMRAI_appu
    SAMRAI_tbox)

target_include_directories( async_br
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/clustering/async_br)

target_compile_definitions(async_br PUBLIC TESTING=1)
