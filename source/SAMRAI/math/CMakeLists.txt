set (math_headers
  ArrayDataBasicOps.h
  ArrayDataMiscellaneousOpsReal.h
  ArrayDataNormOpsComplex.h
  ArrayDataNormOpsInteger.h
  ArrayDataNormOpsReal.h
  HierarchyCellDataOpsComplex.h
  HierarchyCellDataOpsInteger.h
  HierarchyCellDataOpsReal.h
  HierarchyDataOpsComplex.h
  HierarchyDataOpsInteger.h
  HierarchyDataOpsManager.h
  HierarchyDataOpsReal.h
  HierarchyEdgeDataOpsComplex.h
  HierarchyEdgeDataOpsInteger.h
  HierarchyEdgeDataOpsReal.h
  HierarchyFaceDataOpsComplex.h
  HierarchyFaceDataOpsInteger.h
  HierarchyFaceDataOpsReal.h
  HierarchyNodeDataOpsComplex.h
  HierarchyNodeDataOpsInteger.h
  HierarchyNodeDataOpsReal.h
  HierarchySideDataOpsComplex.h
  HierarchySideDataOpsInteger.h
  HierarchySideDataOpsReal.h
  PatchCellDataBasicOps.h
  PatchCellDataMiscellaneousOpsReal.h
  PatchCellDataNormOpsComplex.h
  PatchCellDataNormOpsReal.h
  PatchCellDataOpsComplex.h
  PatchCellDataOpsInteger.h
  PatchCellDataOpsReal.h
  PatchEdgeDataBasicOps.h
  PatchEdgeDataMiscellaneousOpsReal.h
  PatchEdgeDataNormOpsComplex.h
  PatchEdgeDataNormOpsReal.h
  PatchEdgeDataOpsComplex.h
  PatchEdgeDataOpsInteger.h
  PatchEdgeDataOpsReal.h
  PatchFaceDataBasicOps.h
  PatchFaceDataMiscellaneousOpsReal.h
  PatchFaceDataNormOpsComplex.h
  PatchFaceDataNormOpsReal.h
  PatchFaceDataOpsComplex.h
  PatchFaceDataOpsInteger.h
  PatchFaceDataOpsReal.h
  PatchNodeDataBasicOps.h
  PatchNodeDataMiscellaneousOpsReal.h
  PatchNodeDataNormOpsComplex.h
  PatchNodeDataNormOpsReal.h
  PatchNodeDataOpsComplex.h
  PatchNodeDataOpsInteger.h
  PatchNodeDataOpsReal.h
  PatchSideDataBasicOps.h
  PatchSideDataMiscellaneousOpsReal.h
  PatchSideDataNormOpsComplex.h
  PatchSideDataNormOpsReal.h
  PatchSideDataOpsComplex.h
  PatchSideDataOpsInteger.h
  PatchSideDataOpsReal.h)

set (math_sources
  ArrayDataNormOpsComplex.C
  ArrayDataNormOpsInteger.C
  HierarchyCellDataOpsComplex.C
  HierarchyCellDataOpsInteger.C
  HierarchyDataOpsComplex.C
  HierarchyDataOpsInteger.C
  HierarchyDataOpsManager.C
  HierarchyEdgeDataOpsComplex.C
  HierarchyEdgeDataOpsInteger.C
  HierarchyFaceDataOpsComplex.C
  HierarchyFaceDataOpsInteger.C
  HierarchyNodeDataOpsComplex.C
  HierarchyNodeDataOpsInteger.C
  HierarchySideDataOpsComplex.C
  HierarchySideDataOpsInteger.C
  PatchCellDataNormOpsComplex.C
  PatchCellDataOpsComplex.C
  PatchCellDataOpsInteger.C
  PatchEdgeDataNormOpsComplex.C
  PatchEdgeDataOpsComplex.C
  PatchEdgeDataOpsInteger.C
  PatchFaceDataNormOpsComplex.C
  PatchFaceDataOpsComplex.C
  PatchFaceDataOpsInteger.C
  PatchNodeDataNormOpsComplex.C
  PatchNodeDataOpsComplex.C
  PatchNodeDataOpsInteger.C
  PatchSideDataNormOpsComplex.C
  PatchSideDataOpsComplex.C
  PatchSideDataOpsInteger.C)

set (math_depends_on
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_tbox)

blt_add_library(
  NAME SAMRAI_math
  SOURCES ${math_sources}
  HEADERS ${math_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/math
  DEPENDS_ON ${math_depends_on})

target_include_directories( SAMRAI_math
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
