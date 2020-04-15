set (xfer_headers
  BoxGeometryVariableFillPattern.h
  CoarsenAlgorithm.h
  CoarsenClasses.h
  CoarsenCopyTransaction.h
  CoarsenPatchStrategy.h
  CoarsenSchedule.h
  CoarsenTransactionFactory.h
  CompositeBoundaryAlgorithm.h
  CompositeBoundarySchedule.h
  PatchInteriorVariableFillPattern.h
  PatchLevelBorderAndInteriorFillPattern.h
  PatchLevelBorderFillPattern.h
  PatchLevelEnhancedFillPattern.h
  PatchLevelFillPattern.h
  PatchLevelFullFillPattern.h
  PatchLevelInteriorFillPattern.h
  RefineAlgorithm.h
  RefineClasses.h
  RefineCopyTransaction.h
  RefinePatchStrategy.h
  RefineSchedule.h
  RefineScheduleConnectorWidthRequestor.h
  RefineTimeTransaction.h
  RefineTransactionFactory.h
  SingularityPatchStrategy.h
  StandardCoarsenTransactionFactory.h
  StandardRefineTransactionFactory.h
  VariableFillPattern.h)

set (xfer_sources
  BoxGeometryVariableFillPattern.C
  CoarsenAlgorithm.C
  CoarsenClasses.C
  CoarsenCopyTransaction.C
  CoarsenPatchStrategy.C
  CoarsenSchedule.C
  CoarsenTransactionFactory.C
  CompositeBoundaryAlgorithm.C
  CompositeBoundarySchedule.C
  PatchInteriorVariableFillPattern.C
  PatchLevelBorderAndInteriorFillPattern.C
  PatchLevelBorderFillPattern.C
  PatchLevelEnhancedFillPattern.C
  PatchLevelFillPattern.C
  PatchLevelFullFillPattern.C
  PatchLevelInteriorFillPattern.C
  RefineAlgorithm.C
  RefineClasses.C
  RefineCopyTransaction.C
  RefinePatchStrategy.C
  RefineSchedule.C
  RefineScheduleConnectorWidthRequestor.C
  RefineTimeTransaction.C
  RefineTransactionFactory.C
  SingularityPatchStrategy.C
  StandardCoarsenTransactionFactory.C
  StandardRefineTransactionFactory.C
  VariableFillPattern.C)

set (xfer_depends 
   SAMRAI_hier
   SAMRAI_tbox)

if (ENABLE_CUDA)
  set(xfer_depends ${xfer_depends} cuda)
endif ()

blt_add_library(
  NAME SAMRAI_xfer
  SOURCES ${xfer_sources}
  HEADERS ${xfer_headers}
  DEPENDS_ON ${xfer_depends})

target_include_directories(
  SAMRAI_xfer
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SAMRAI_xfer
  EXPORT SAMRAITargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${xfer_headers}
  DESTINATION include/SAMRAI/xfer)
