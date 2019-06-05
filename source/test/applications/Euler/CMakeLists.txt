set (euler_sources
  Euler.C
  main.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/avgint2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/avgint3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flaten1d2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flaten1d3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flux2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flux3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/fluxcorner3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/grad2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/grad3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stable2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stable3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stufprobc.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace1d2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace1d3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace3d.f)

process_m4(NAME fortran/avgint2d)
process_m4(NAME fortran/avgint3d)
process_m4(NAME fortran/flaten1d2d)
process_m4(NAME fortran/flaten1d3d)
process_m4(NAME fortran/flux2d)
process_m4(NAME fortran/flux3d)
process_m4(NAME fortran/fluxcorner3d)
process_m4(NAME fortran/grad2d)
process_m4(NAME fortran/grad3d)
process_m4(NAME fortran/init2d)
process_m4(NAME fortran/init3d)
process_m4(NAME fortran/stable2d)
process_m4(NAME fortran/stable3d)
process_m4(NAME fortran/stufprobc)
process_m4(NAME fortran/trace1d2d)
process_m4(NAME fortran/trace1d3d)
process_m4(NAME fortran/trace2d)
process_m4(NAME fortran/trace3d)

blt_add_executable(
  NAME euler
  SOURCES ${euler_sources}
  DEPENDS_ON
    SAMRAI_algs
    SAMRAI_appu
    SAMRAI_geom
    SAMRAI_hier
    SAMRAI_math
    SAMRAI_mesh
    SAMRAI_pdat
    SAMRAI_tbox
    SAMRAI_xfer)

target_include_directories( euler
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/applications/Euler)

set_target_properties(euler PROPERTIES LINKER_LANGUAGE CXX)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME euler
  EXECUTABLE euler
  INPUTS ${test_inputs}
  PARALLEL TRUE)

