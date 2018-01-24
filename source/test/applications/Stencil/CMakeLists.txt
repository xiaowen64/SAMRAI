set (stencil_sources
  Stencil.C
  main.C
)

blt_add_executable(
  NAME stencil
  SOURCES ${stencil_sources}
  DEPENDS_ON
    raja
    cuda
    openmp
    SAMRAI_algs
    SAMRAI_appu
    SAMRAI_geom
    SAMRAI_hier
    SAMRAI_math
    SAMRAI_mesh
    SAMRAI_pdat
    SAMRAI_tbox
    SAMRAI_xfer)

target_include_directories( stencil
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/applications/Stencil)

set_target_properties(stencil PROPERTIES LINKER_LANGUAGE CXX)
