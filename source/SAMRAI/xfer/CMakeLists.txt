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

set ( xfer_depends_on
  SAMRAI_hier)

blt_add_library(
  NAME SAMRAI_xfer
  SOURCES ${xfer_sources}
  HEADERS ${xfer_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/xfer
  DEPENDS_ON ${xfer_depends_on})

target_include_directories(
  SAMRAI_xfer
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
