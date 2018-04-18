set (solv_headers
  CartesianRobinBcHelper.h
  CellPoissonFACOps.h
  CellPoissonFACSolver.h
  CellPoissonHypreSolver.h
  CVODEAbstractFunctions.h
  CVODESolver.h
  FACOperatorStrategy.h
  FACPreconditioner.h
  GhostCellRobinBcCoefs.h
  KINSOL_SAMRAIContext.h
  KINSOLAbstractFunctions.h
  KINSOLSolver.h
  LocationIndexRobinBcCoefs.h
  NonlinearSolverStrategy.h
  nvector.h
  PETSc_SAMRAIVectorReal.h
  PETSc_SAMRAIVectorReal.C
  PETScAbstractVectorReal.h
  PETScAbstractVectorReal.C
  PoissonSpecifications.h
  RobinBcCoefStrategy.h
  SAMRAIVectorReal.h
  SAMRAIVectorReal.C
  SimpleCellRobinBcCoefs.h
  SNES_SAMRAIContext.h
  SNESAbstractFunctions.h
  solv_NVector.h
  Sundials_SAMRAIVector.h
  SundialsAbstractVector.h)

set_source_files_properties(
  PETScAbstractVectorReal.C
  PETSc_SAMRAIVectorReal.C
  SAMRAIVectorReal.C
  PROPERTIES HEADER_FILE_ONLY TRUE)

set (solv_sources
  CVODEAbstractFunctions.C
  CVODESolver.C
  CartesianRobinBcHelper.C
  CellPoissonFACOps.C
  CellPoissonFACSolver.C
  CellPoissonHypreSolver.C
  FACOperatorStrategy.C
  FACPreconditioner.C
  GhostCellRobinBcCoefs.C
  KINSOLAbstractFunctions.C
  KINSOLSolver.C
  KINSOL_SAMRAIContext.C
  LocationIndexRobinBcCoefs.C
  NonlinearSolverStrategy.C
  PoissonSpecifications.C
  RobinBcCoefStrategy.C
  SNESAbstractFunctions.C
  SNES_SAMRAIContext.C
  SimpleCellRobinBcCoefs.C
  SundialsAbstractVector.C
  Sundials_SAMRAIVector.C
  solv_NVector.C
  fortran/solv_cartesianrobinbchelper2d.f
  fortran/solv_cartesianrobinbchelper3d.f
  fortran/solv_cellpoissonfacops2d.f
  fortran/solv_cellpoissonfacops3d.f
  fortran/solv_cellpoissonhypresolver2d.f
  fortran/solv_cellpoissonhypresolver3d.f)

set_source_files_properties(
  fortran/solv_cartesianrobinbchelper2d.f
  fortran/solv_cartesianrobinbchelper3d.f
  fortran/solv_cellpoissonfacops2d.f
  fortran/solv_cellpoissonfacops3d.f
  fortran/solv_cellpoissonhypresolver2d.f
  fortran/solv_cellpoissonhypresolver3d.f
  PROPERTIES
  LANGUAGE
  Fortran)

set (solv_depends_on
  SAMRAI_hier
  SAMRAI_math
  SAMRAI_tbox)

if (ENABLE_CUDA)
  set (solv_depends_on
    ${solv_depends_on}
    cuda)
endif ()

blt_add_library(
  NAME SAMRAI_solv
  SOURCES ${solv_sources}
  HEADERS ${solv_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/solv
  DEPENDS_ON ${solv_depends_on})

install(TARGETS
  SAMRAI_solv
  RUNTIME DESTINATION lib
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

set_target_properties(SAMRAI_solv PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(
  SAMRAI_solv
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
