set (multiblock_sources
  main.C
  MblkGeometry.C
  MblkLinAdv.C
  SkeletonOutersideDoubleWeightedAverage.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stufprobc.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/trace3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/grad2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/grad3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flux2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/flux3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/fluxcorner3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stable2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/stable3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/init3d.f)

process_m4(NAME fortran/stufprobc)
process_m4(NAME fortran/trace1d)
process_m4(NAME fortran/trace2d)
process_m4(NAME fortran/trace3d)
process_m4(NAME fortran/grad2d)
process_m4(NAME fortran/grad3d)
process_m4(NAME fortran/flux2d)
process_m4(NAME fortran/flux3d)
process_m4(NAME fortran/fluxcorner3d)
process_m4(NAME fortran/stable2d)
process_m4(NAME fortran/stable3d)
process_m4(NAME fortran/init2d)
process_m4(NAME fortran/init3d)

blt_add_executable(
  NAME multiblock
  SOURCES ${multiblock_sources}
  DEPENDS_ON
    testlib
    SAMRAI)

  #target_compile_definitions(multiblock PUBLIC TESTING=1)
