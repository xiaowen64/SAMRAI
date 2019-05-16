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

set(convdiff_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(convdiff_depends ${convdiff_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(convdiff_depends ${convdiff_depends} cuda)
endif ()

blt_add_executable(
  NAME convdiff
  SOURCES ${convdiff_sources}
  DEPENDS_ON ${convdiff_depends})

target_include_directories(
  convdiff
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/applications/ConvDiff)

set_target_properties(convdiff PROPERTIES LINKER_LANGUAGE CXX)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME convdiff
  EXECUTABLE convdiff
  INPUTS ${test_inputs}
  PARALLEL TRUE)

