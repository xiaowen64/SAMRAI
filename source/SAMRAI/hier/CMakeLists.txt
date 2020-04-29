set (hier_headers
  AssumedPartition.h
  AssumedPartitionBox.h
  BaseConnectorAlgorithm.h
  BaseGridGeometry.h
  BlockId.h
  BlueprintUtils.h
  BlueprintUtilsStrategy.h
  BoundaryBox.h
  BoundaryBoxUtils.h
  BoundaryLookupTable.h
  Box.h
  BoxContainer.h
  BoxContainerSingleBlockIterator.h
  BoxContainerSingleOwnerIterator.h
  BoxContainerUtils.h
  BoxGeometry.h
  BoxId.h
  BoxLevel.h
  BoxLevelConnectorUtils.h
  BoxLevelHandle.h
  BoxLevelStatistics.h
  BoxNeighborhoodCollection.h
  BoxOverlap.h
  BoxTree.h
  BoxUtilities.h
  CoarseFineBoundary.h
  CoarsenOperator.h
  ComponentSelector.h
  Connector.h
  ConnectorStatistics.h
  FlattenedHierarchy.h
  GlobalId.h
  HierarchyNeighbors.h
  Index.h
  IntVector.h
  LocalId.h
  MappingConnector.h
  MappingConnectorAlgorithm.h
  MultiblockBoxTree.h
  OverlapConnectorAlgorithm.h
  Patch.h
  PatchBoundaries.h
  PatchData.h
  PatchDataFactory.h
  PatchDataRestartManager.h
  PatchDescriptor.h
  PatchFactory.h
  PatchGeometry.h
  PatchHierarchy.h
  PatchLevel.h
  PatchLevelFactory.h
  PeriodicId.h
  PeriodicShiftCatalog.h
  PersistentOverlapConnectors.h
  ProcessorMapping.h
  RealBoxConstIterator.h
  RefineOperator.h
  SequentialLocalIdGenerator.h
  SingularityFinder.h
  TimeInterpolateOperator.h
  TransferOperatorRegistry.h
  Transformation.h
  UncoveredBoxIterator.h
  Variable.h
  VariableContext.h
  VariableDatabase.h)

set (hier_sources
  AssumedPartition.C
  AssumedPartitionBox.C
  BaseConnectorAlgorithm.C
  BaseGridGeometry.C
  BlockId.C
  BlueprintUtils.C
  BlueprintUtilsStrategy.C
  BoundaryBox.C
  BoundaryBoxUtils.C
  BoundaryLookupTable.C
  Box.C
  BoxContainer.C
  BoxContainerSingleBlockIterator.C
  BoxContainerSingleOwnerIterator.C
  BoxContainerUtils.C
  BoxGeometry.C
  BoxId.C
  BoxLevel.C
  BoxLevelConnectorUtils.C
  BoxLevelHandle.C
  BoxLevelStatistics.C
  BoxNeighborhoodCollection.C
  BoxOverlap.C
  BoxTree.C
  BoxUtilities.C
  CoarseFineBoundary.C
  CoarsenOperator.C
  ComponentSelector.C
  Connector.C
  ConnectorStatistics.C
  FlattenedHierarchy.C
  GlobalId.C
  HierarchyNeighbors.C
  Index.C
  IntVector.C
  LocalId.C
  MappingConnector.C
  MappingConnectorAlgorithm.C
  MultiblockBoxTree.C
  OverlapConnectorAlgorithm.C
  Patch.C
  PatchBoundaries.C
  PatchData.C
  PatchDataFactory.C
  PatchDataRestartManager.C
  PatchDescriptor.C
  PatchFactory.C
  PatchGeometry.C
  PatchHierarchy.C
  PatchLevel.C
  PatchLevelFactory.C
  PeriodicId.C
  PeriodicShiftCatalog.C
  PersistentOverlapConnectors.C
  ProcessorMapping.C
  RealBoxConstIterator.C
  RefineOperator.C
  SingularityFinder.C
  TimeInterpolateOperator.C
  TransferOperatorRegistry.C
  Transformation.C
  UncoveredBoxIterator.C
  Variable.C
  VariableContext.C
  VariableDatabase.C)

set (hier_depends SAMRAI_tbox)

if (ENABLE_MPI)
  set (hier_depends ${hier_depends} mpi)
endif ()

blt_add_library(
  NAME SAMRAI_hier
  SOURCES ${hier_sources}
  HEADERS ${hier_headers}
  DEPENDS_ON ${hier_depends})

target_include_directories(
  SAMRAI_hier
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SAMRAI_hier
  EXPORT SAMRAITargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${hier_headers}
  DESTINATION include/SAMRAI/hier)
