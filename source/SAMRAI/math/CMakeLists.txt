set (math_headers
  ArrayDataBasicOps.C
  ArrayDataBasicOps.h
  ArrayDataMiscellaneousOpsReal.C
  ArrayDataMiscellaneousOpsReal.h
  ArrayDataNormOpsComplex.h
  ArrayDataNormOpsInteger.h
  ArrayDataNormOpsReal.C
  ArrayDataNormOpsReal.h
  HierarchyCellDataOpsComplex.h
  HierarchyCellDataOpsInteger.h
  HierarchyCellDataOpsReal.C
  HierarchyCellDataOpsReal.h
  HierarchyDataOpsComplex.h
  HierarchyDataOpsInteger.h
  HierarchyDataOpsManager.h
  HierarchyDataOpsReal.C
  HierarchyDataOpsReal.h
  HierarchyEdgeDataOpsComplex.h
  HierarchyEdgeDataOpsInteger.h
  HierarchyEdgeDataOpsReal.C
  HierarchyEdgeDataOpsReal.h
  HierarchyFaceDataOpsComplex.h
  HierarchyFaceDataOpsInteger.h
  HierarchyFaceDataOpsReal.C
  HierarchyFaceDataOpsReal.h
  HierarchyNodeDataOpsComplex.h
  HierarchyNodeDataOpsInteger.h
  HierarchyNodeDataOpsReal.C
  HierarchyNodeDataOpsReal.h
  HierarchySideDataOpsComplex.h
  HierarchySideDataOpsInteger.h
  HierarchySideDataOpsReal.C
  HierarchySideDataOpsReal.h
  PatchCellDataBasicOps.C
  PatchCellDataBasicOps.h
  PatchCellDataMiscellaneousOpsReal.C
  PatchCellDataMiscellaneousOpsReal.h
  PatchCellDataNormOpsComplex.h
  PatchCellDataNormOpsReal.C
  PatchCellDataNormOpsReal.h
  PatchCellDataOpsComplex.h
  PatchCellDataOpsInteger.h
  PatchCellDataOpsReal.C
  PatchCellDataOpsReal.h
  PatchEdgeDataBasicOps.C
  PatchEdgeDataBasicOps.h
  PatchEdgeDataMiscellaneousOpsReal.C
  PatchEdgeDataMiscellaneousOpsReal.h
  PatchEdgeDataNormOpsComplex.h
  PatchEdgeDataNormOpsReal.C
  PatchEdgeDataNormOpsReal.h
  PatchEdgeDataOpsComplex.h
  PatchEdgeDataOpsInteger.h
  PatchEdgeDataOpsReal.C
  PatchEdgeDataOpsReal.h
  PatchFaceDataBasicOps.C
  PatchFaceDataBasicOps.h
  PatchFaceDataMiscellaneousOpsReal.C
  PatchFaceDataMiscellaneousOpsReal.h
  PatchFaceDataNormOpsComplex.h
  PatchFaceDataNormOpsReal.C
  PatchFaceDataNormOpsReal.h
  PatchFaceDataOpsComplex.h
  PatchFaceDataOpsInteger.h
  PatchFaceDataOpsReal.C
  PatchFaceDataOpsReal.h
  PatchNodeDataBasicOps.C
  PatchNodeDataBasicOps.h
  PatchNodeDataMiscellaneousOpsReal.C
  PatchNodeDataMiscellaneousOpsReal.h
  PatchNodeDataNormOpsComplex.h
  PatchNodeDataNormOpsReal.C
  PatchNodeDataNormOpsReal.h
  PatchNodeDataOpsComplex.h
  PatchNodeDataOpsInteger.h
  PatchNodeDataOpsReal.C
  PatchNodeDataOpsReal.h
  PatchSideDataBasicOps.C
  PatchSideDataBasicOps.h
  PatchSideDataMiscellaneousOpsReal.C
  PatchSideDataMiscellaneousOpsReal.h
  PatchSideDataNormOpsComplex.h
  PatchSideDataNormOpsReal.C
  PatchSideDataNormOpsReal.h
  PatchSideDataOpsComplex.h
  PatchSideDataOpsInteger.h
  PatchSideDataOpsReal.C
  PatchSideDataOpsReal.h)

set_source_files_properties(
  ArrayDataBasicOps.C
  ArrayDataMiscellaneousOpsReal.C
  ArrayDataNormOpsReal.C
  HierarchyCellDataOpsReal.C
  HierarchyDataOpsReal.C
  HierarchyEdgeDataOpsReal.C
  HierarchyFaceDataOpsReal.C
  HierarchyNodeDataOpsReal.C
  HierarchySideDataOpsReal.C
  PatchCellDataBasicOps.C
  PatchCellDataMiscellaneousOpsReal.C
  PatchCellDataNormOpsReal.C
  PatchCellDataOpsReal.C
  PatchEdgeDataBasicOps.C
  PatchEdgeDataMiscellaneousOpsReal.C
  PatchEdgeDataNormOpsReal.C
  PatchEdgeDataOpsReal.C
  PatchFaceDataBasicOps.C
  PatchFaceDataMiscellaneousOpsReal.C
  PatchFaceDataNormOpsReal.C
  PatchFaceDataOpsReal.C
  PatchNodeDataBasicOps.C
  PatchNodeDataMiscellaneousOpsReal.C
  PatchNodeDataNormOpsReal.C
  PatchNodeDataOpsReal.C
  PatchSideDataBasicOps.C
  PatchSideDataMiscellaneousOpsReal.C
  PatchSideDataNormOpsReal.C
  PatchSideDataOpsReal.C
  PROPERTIES HEADER_FILE_ONLY True )

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

set (math_depends
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_tbox)

if (ENABLE_CUDA)
  set (math_depends ${math_depends} cuda)
endif ()

blt_add_library(
  NAME SAMRAI_math
  SOURCES ${math_sources}
  HEADERS ${math_headers}
  DEPENDS_ON ${math_depends})

target_include_directories( SAMRAI_math
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SAMRAI_math
  EXPORT SAMRAITargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${math_headers}
  DESTINATION include/SAMRAI/math)
