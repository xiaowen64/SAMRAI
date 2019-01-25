set (facstaticrefinement_sources
  FACPoisson.C
  main.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/facpoisson2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/facpoisson3d.f)

process_m4(NAME fortran/facpoisson2d)
process_m4(NAME fortran/facpoisson3d)

set (facstaticrefinement_depends_on
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_mesh
  SAMRAI_solv
  SAMRAI_appu
  SAMRAI_xfer
  SAMRAI_geom
  SAMRAI_tbox)

blt_add_executable(
  NAME facstaticrefinement
  SOURCES ${facstaticrefinement_sources}
  DEPENDS_ON
    ${facstaticrefinement_depends_on})

target_compile_definitions(facstaticrefinement PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME facstaticrefinement
  EXECUTABLE facstaticrefinement
  INPUTS ${test_inputs}
  PARALLEL TRUE)
