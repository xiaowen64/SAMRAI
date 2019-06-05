set (testlib_headers
  AutoTester.h
  DerivedVisOwnerData.h
  GeomUtilsAMR.h
  get-input-filename.h
  MblkHyperbolicLevelIntegrator.h
  MblkHyperbolicPatchStrategy.h
  MeshGenerationStrategy.h
  ShrunkenLevelGenerator.h
  SinusoidalFrontGenerator.h
  SkeletonBoundaryUtilities2.h
  SkeletonBoundaryUtilities3.h
  SphericalShellGenerator.h)

set (testlib_sources
  AutoTester.C
  DerivedVisOwnerData.C
  MblkHyperbolicLevelIntegrator.C
  MblkHyperbolicPatchStrategy.C
  SkeletonBoundaryUtilities2.C
  SkeletonBoundaryUtilities3.C
  SkeletonCellDoubleConservativeLinearRefine.C
  SkeletonCellDoubleWeightedAverage.C
  MeshGenerationStrategy.C
  SinusoidalFrontGenerator.C
  SphericalShellGenerator.C
  ShrunkenLevelGenerator.C
  get-input-filename.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/skelbdry2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/skelbdry3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/upfluxsum2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/upfluxsum3d.f)

set (testlib_depends_on
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_tbox
  )

process_m4(NAME fortran/skelbdry2d)
process_m4(NAME fortran/skelbdry3d)
process_m4(NAME fortran/upfluxsum2d)
process_m4(NAME fortran/upfluxsum3d)

blt_add_library( NAME testlib
  SOURCES ${testlib_sources}
  HEADERS ${testlib_headers}
  DEPENDS_ON ${testlib_depends_on})

set_target_properties(testlib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories( testlib
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/testlib)

target_compile_definitions(testlib PRIVATE TESTING=1)
