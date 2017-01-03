set (pdat_headers
  ArrayData.h
  ArrayData.C
  ArrayDataAccess.h
  ArrayDataIterator.h
  ArrayDataOperationUtilities.h
  ArrayDataOperationUtilities.C
  CellComplexConstantRefine.h
  CellComplexLinearTimeInterpolateOp.h
  CellData.h
  CellData.C
  CellDataFactory.h
  CellDataFactory.C
  CellDoubleConstantRefine.h
  CellDoubleLinearTimeInterpolateOp.h
  CellFloatConstantRefine.h
  CellFloatLinearTimeInterpolateOp.h
  CellGeometry.h
  CellIndex.h
  CellIntegerConstantRefine.h
  CellIterator.h
  CellOverlap.h
  CellVariable.h
  CellVariable.C
  CopyOperation.h
  CopyOperation.C
  DoubleAttributeId.h
  EdgeComplexConstantRefine.h
  EdgeComplexLinearTimeInterpolateOp.h
  EdgeData.h
  EdgeDataFactory.h
  EdgeDoubleConstantRefine.h
  EdgeDoubleLinearTimeInterpolateOp.h
  EdgeFloatConstantRefine.h
  EdgeFloatLinearTimeInterpolateOp.h
  EdgeGeometry.h
  EdgeIndex.h
  EdgeIntegerConstantRefine.h
  EdgeIterator.h
  EdgeOverlap.h
  EdgeVariable.h
  FaceComplexConstantRefine.h
  FaceComplexLinearTimeInterpolateOp.h
  FaceData.h
  FaceDataFactory.h
  FaceDoubleConstantRefine.h
  FaceDoubleLinearTimeInterpolateOp.h
  FaceFloatConstantRefine.h
  FaceFloatLinearTimeInterpolateOp.h
  FaceGeometry.h
  FaceIndex.h
  FaceIntegerConstantRefine.h
  FaceIterator.h
  FaceOverlap.h
  FaceVariable.h
  FirstLayerCellNoCornersVariableFillPattern.h
  FirstLayerCellVariableFillPattern.h
  FirstLayerEdgeVariableFillPattern.h
  FirstLayerNodeVariableFillPattern.h
  FirstLayerSideVariableFillPattern.h
  IndexData.h
  IndexDataFactory.h
  IndexVariable.h
  IntegerAttributeId.h
  MDA_Access.h
  NodeComplexInjection.h
  NodeComplexLinearTimeInterpolateOp.h
  NodeData.h
  NodeData.C
  NodeDataFactory.h
  NodeDoubleInjection.h
  NodeDoubleLinearTimeInterpolateOp.h
  NodeFloatInjection.h
  NodeFloatLinearTimeInterpolateOp.h
  NodeGeometry.h
  NodeIndex.h
  NodeIntegerInjection.h
  NodeIterator.h
  NodeOverlap.h
  NodeVariable.h
  NodeVariable.C
  OuteredgeData.h
  OuteredgeDataFactory.h
  OuteredgeGeometry.h
  OuteredgeVariable.h
  OuterfaceComplexConstantRefine.h
  OuterfaceComplexLinearTimeInterpolateOp.h
  OuterfaceData.h
  OuterfaceDataFactory.h
  OuterfaceDoubleConstantRefine.h
  OuterfaceDoubleLinearTimeInterpolateOp.h
  OuterfaceFloatConstantRefine.h
  OuterfaceFloatLinearTimeInterpolateOp.h
  OuterfaceGeometry.h
  OuterfaceIntegerConstantRefine.h
  OuterfaceVariable.h
  OuternodeData.h
  OuternodeDataFactory.h
  OuternodeDoubleInjection.h
  OuternodeGeometry.h
  OuternodeVariable.h
  OutersideComplexLinearTimeInterpolateOp.h
  OutersideData.h
  OutersideDataFactory.h
  OutersideDoubleLinearTimeInterpolateOp.h
  OutersideFloatLinearTimeInterpolateOp.h
  OutersideGeometry.h
  OutersideVariable.h
  SecondLayerNodeNoCornersVariableFillPattern.h
  SecondLayerNodeVariableFillPattern.h
  SideComplexConstantRefine.h
  SideComplexLinearTimeInterpolateOp.h
  SideData.C
  SideData.h
  SideDataFactory.h
  SideDoubleConstantRefine.h
  SideDoubleLinearTimeInterpolateOp.h
  SideFloatConstantRefine.h
  SideFloatLinearTimeInterpolateOp.h
  SideGeometry.h
  SideIndex.h
  SideIntegerConstantRefine.h
  SideIterator.h
  SideOverlap.h
  SideVariable.C
  SideVariable.h
  SparseData.h
  SparseDataFactory.h
  SparseDataVariable.h
  SumOperation.h
  SumOperation.C)

set_source_files_properties(
  ArrayData.C
  ArrayDataOperationUtilities.C
  CellData.C
  CellDataFactory.C
  CellVariable.C
  CopyOperation.C
  NodeData.C
  NodeVariable.C
  SideData.C
  SideVariable.C
  SumOperation.C
  PROPERTIES HEADER_FILE_ONLY TRUE)


set (pdat_sources
  ArrayDataIterator.C
  CellComplexConstantRefine.C
  CellComplexLinearTimeInterpolateOp.C
  CellDoubleConstantRefine.C
  CellDoubleLinearTimeInterpolateOp.C
  CellFloatConstantRefine.C
  CellFloatLinearTimeInterpolateOp.C
  CellGeometry.C
  CellIndex.C
  CellIntegerConstantRefine.C
  CellIterator.C
  CellOverlap.C
  DoubleAttributeId.C
  EdgeComplexConstantRefine.C
  EdgeComplexLinearTimeInterpolateOp.C
  EdgeDoubleConstantRefine.C
  EdgeDoubleLinearTimeInterpolateOp.C
  EdgeFloatConstantRefine.C
  EdgeFloatLinearTimeInterpolateOp.C
  EdgeGeometry.C
  EdgeIndex.C
  EdgeIntegerConstantRefine.C
  EdgeIterator.C
  EdgeOverlap.C
  FaceComplexConstantRefine.C
  FaceComplexLinearTimeInterpolateOp.C
  FaceDoubleConstantRefine.C
  FaceDoubleLinearTimeInterpolateOp.C
  FaceFloatConstantRefine.C
  FaceFloatLinearTimeInterpolateOp.C
  FaceGeometry.C
  FaceIndex.C
  FaceIntegerConstantRefine.C
  FaceIterator.C
  FaceOverlap.C
  FirstLayerCellNoCornersVariableFillPattern.C
  FirstLayerCellVariableFillPattern.C
  FirstLayerEdgeVariableFillPattern.C
  FirstLayerNodeVariableFillPattern.C
  FirstLayerSideVariableFillPattern.C
  IntegerAttributeId.C
  NodeComplexInjection.C
  NodeComplexLinearTimeInterpolateOp.C
  NodeDoubleInjection.C
  NodeDoubleLinearTimeInterpolateOp.C
  NodeFloatInjection.C
  NodeFloatLinearTimeInterpolateOp.C
  NodeGeometry.C
  NodeIndex.C
  NodeIntegerInjection.C
  NodeIterator.C
  NodeOverlap.C
  OuteredgeGeometry.C
  OuterfaceComplexConstantRefine.C
  OuterfaceComplexLinearTimeInterpolateOp.C
  OuterfaceDoubleConstantRefine.C
  OuterfaceDoubleLinearTimeInterpolateOp.C
  OuterfaceFloatConstantRefine.C
  OuterfaceFloatLinearTimeInterpolateOp.C
  OuterfaceGeometry.C
  OuterfaceIntegerConstantRefine.C
  OuternodeDoubleInjection.C
  OuternodeGeometry.C
  OutersideComplexLinearTimeInterpolateOp.C
  OutersideDoubleLinearTimeInterpolateOp.C
  OutersideFloatLinearTimeInterpolateOp.C
  OutersideGeometry.C
  SecondLayerNodeNoCornersVariableFillPattern.C
  SecondLayerNodeVariableFillPattern.C
  SideComplexConstantRefine.C
  SideComplexLinearTimeInterpolateOp.C
  SideDoubleConstantRefine.C
  SideDoubleLinearTimeInterpolateOp.C
  SideFloatConstantRefine.C
  SideFloatLinearTimeInterpolateOp.C
  SideGeometry.C
  SideIndex.C
  SideIntegerConstantRefine.C
  SideIterator.C
  SideOverlap.C
  SparseData.C
  SparseDataFactory.C
  SparseDataVariable.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_concoarsen1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_concoarsen2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_concoarsen3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_conrefine1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_conrefine2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_conrefine3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_dbugfort1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_dbugfort2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_dbugfort3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_lintimint1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_lintimint2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/pdat_lintimint3d.f)

set (pdat_depends_on
  SAMRAI_hier
  SAMRAI_tbox)

process_m4(NAME fortran/pdat_concoarsen1d)
process_m4(NAME fortran/pdat_concoarsen2d)
process_m4(NAME fortran/pdat_concoarsen3d)
process_m4(NAME fortran/pdat_conrefine1d)
process_m4(NAME fortran/pdat_conrefine2d)
process_m4(NAME fortran/pdat_conrefine3d)
process_m4(NAME fortran/pdat_dbugfort1d)
process_m4(NAME fortran/pdat_dbugfort2d)
process_m4(NAME fortran/pdat_dbugfort3d)
process_m4(NAME fortran/pdat_lintimint1d)
process_m4(NAME fortran/pdat_lintimint2d)
process_m4(NAME fortran/pdat_lintimint3d)

blt_add_library(
  NAME SAMRAI_pdat
  SOURCES ${pdat_sources}
  HEADERS ${pdat_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/pdat
  DEPENDS_ON ${pdat_depends_on})

set_target_properties(SAMRAI_pdat PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(
  SAMRAI_pdat
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
