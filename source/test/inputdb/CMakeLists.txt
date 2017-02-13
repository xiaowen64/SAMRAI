set (inputdb_sources
        inputdb.C)

blt_add_executable(
        NAME inputdb
        SOURCES ${inputdb_sources}
        DEPENDS_ON
          SAMRAI_tbox)

target_compile_definitions(
        inputdb
        PRIVATE
        NDIM=2 TESTING=1)