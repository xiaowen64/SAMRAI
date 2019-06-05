set (treesearch_sources
  main.C)

blt_add_executable(
  NAME treesearch
  SOURCES ${treesearch_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_tbox)

target_compile_definitions(treesearch PUBLIC TESTING=1)
