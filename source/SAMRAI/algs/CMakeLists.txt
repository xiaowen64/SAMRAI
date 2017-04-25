set (algs_headers
  HyperbolicLevelIntegrator.h
  HyperbolicPatchStrategy.h
  ImplicitEquationStrategy.h
  ImplicitIntegrator.h
  MethodOfLinesIntegrator.h
  MethodOfLinesPatchStrategy.h
  OuteredgeSumTransaction.h
  OuteredgeSumTransactionFactory.h
  OuternodeSumTransaction.h
  OuternodeSumTransactionFactory.h
  PatchBoundaryEdgeSum.h
  PatchBoundaryNodeSum.h
  TimeRefinementIntegrator.h
  TimeRefinementIntegratorConnectorWidthRequestor.h
  TimeRefinementLevelStrategy.h)

set (algs_sources
  HyperbolicLevelIntegrator.C
  HyperbolicPatchStrategy.C
  ImplicitEquationStrategy.C
  ImplicitIntegrator.C
  MethodOfLinesIntegrator.C
  MethodOfLinesPatchStrategy.C
  OuteredgeSumTransaction.C
  OuteredgeSumTransactionFactory.C
  OuternodeSumTransaction.C
  OuternodeSumTransactionFactory.C
  PatchBoundaryEdgeSum.C
  PatchBoundaryNodeSum.C
  TimeRefinementIntegrator.C
  TimeRefinementIntegratorConnectorWidthRequestor.C
  TimeRefinementLevelStrategy.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/algs_nodeouternodeops2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/algs_nodeouternodeops3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/algs_upfluxsum1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/algs_upfluxsum2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/algs_upfluxsum3d.f)

set (algs_depends_on
  SAMRAI_hier
  SAMRAI_tbox
  SAMRAI_pdat
  boost)

if (ENABLE_CUDA)
  set (algs_depends_on
    ${algs_depends_on}
    cuda)
endif ()

process_m4(NAME fortran/algs_nodeouternodeops2d)
process_m4(NAME fortran/algs_nodeouternodeops3d)
process_m4(NAME fortran/algs_upfluxsum1d)
process_m4(NAME fortran/algs_upfluxsum2d)
process_m4(NAME fortran/algs_upfluxsum3d)

blt_add_library( NAME SAMRAI_algs
  SOURCES ${algs_sources}
  HEADERS ${algs_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/algs
  DEPENDS_ON ${algs_depends_on})

set_target_properties(SAMRAI_algs PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories( SAMRAI_algs
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
