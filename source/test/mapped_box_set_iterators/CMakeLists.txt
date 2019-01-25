set (mapped_box_sources
        main-iterators.C)

blt_add_executable(
        NAME mapped_box_set_iterators
        SOURCES ${mapped_box_sources}
        DEPENDS_ON
          SAMRAI_hier
          SAMRAI_tbox)

blt_add_test(
  NAME mapped_box_set_iterators
  COMMAND mapped_box_set_iterators)
