set (hypre_sources
  main.C
  HyprePoisson.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/hyprepoisson2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/hyprepoisson3d.f)

process_m4(NAME fortran/hyprepoisson2d)
process_m4(NAME fortran/hyprepoisson3d)

set (hypre_depends_on
    SAMRAI_solv
    SAMRAI_pdat
    SAMRAI_mesh
    SAMRAI_geom
    SAMRAI_tbox
    SAMRAI_appu)

if (HAVE_HYPRE)
  set (hypre_depends_on
    ${hypre_depends_on}; hypre)
endif ()

blt_add_executable(
  NAME hypre-test
  SOURCES ${hypre_sources}
  DEPENDS_ON
    ${hypre_depends_on})

target_compile_definitions(hypre-test PUBLIC TESTING=1)
