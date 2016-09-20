set (convdiff_sources
  ConvDiff.C
  main.C
  MainRestartData.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/compute_rhs2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/compute_rhs3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/rkstep2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/rkstep3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/tag_cells2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/tag_cells3d.f)

process_m4(NAME fortran/compute_rhs2d)
process_m4(NAME fortran/compute_rhs3d)
process_m4(NAME fortran/init2d)
process_m4(NAME fortran/init3d)
process_m4(NAME fortran/rkstep2d)
process_m4(NAME fortran/rkstep3d)
process_m4(NAME fortran/tag_cells2d)
process_m4(NAME fortran/tag_cells3d)

blt_add_executable(
  NAME convdiff
  SOURCES ${convdiff_sources}
  DEPENDS_ON
  hdf5
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_geom
  SAMRAI_hier
  SAMRAI_math
  SAMRAI_mesh
  SAMRAI_pdat
  SAMRAI_tbox
  SAMRAI_xfer)

target_include_directories(
  convdiff
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/applications/ConvDiff)

set_target_properties(convdiff PROPERTIES LINKER_LANGUAGE CXX)
