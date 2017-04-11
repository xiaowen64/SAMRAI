set (patchbdrysum_sources
  main.C
  HierSumTest.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/edgeops2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/edgeops3d.f)

process_m4(NAME fortran/edgeops2d)
process_m4(NAME fortran/edgeops3d)

set (patchbdrysum_depends_on
    SAMRAI_hier
    SAMRAI_mesh
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_algs
    SAMRAI_appu
    SAMRAI_tbox)

blt_add_executable(
  NAME patchbdrysum-test
  SOURCES ${patchbdrysum_sources}
  DEPENDS_ON
    ${patchbdrysum_depends_on})

target_compile_definitions(patchbdrysum-test PUBLIC TESTING=1)
