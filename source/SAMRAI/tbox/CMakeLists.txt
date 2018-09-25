set ( tbox_headers
  Array.h
  Array.C
  AsyncCommGroup.h
  AsyncCommPeer.h
  AsyncCommPeer.C
  AsyncCommStage.h
  BalancedDepthFirstTree.h
  BreadthFirstRankTree.h
  CenteredRankTree.h
  Clock.h
  CommGraphWriter.h
  Complex.h
  Database.h
  DatabaseBox.h
  DatabaseFactory.h
  Dimension.h
  Grammar.h
  HDFDatabase.h
  HDFDatabaseFactory.h
  IEEE.h
  InputDatabase.h
  InputManager.h
  IOStream.h
  Logger.h
  MathUtilities.h
  MathUtilities.C
  MemoryDatabase.h
  MemoryDatabaseFactory.h
  MemoryUtilities.h
  MessageStream.h
  NullDatabase.h
  OpenMPUtilities.h
  ParallelBuffer.h
  Parser.h
  PIO.h
  RankGroup.h
  RankTreeStrategy.h
  ReferenceCounter.h
  RestartManager.h
  SAMRAI_MPI.h
  SAMRAIManager.h
  Schedule.h
  Serializable.h
  SiloDatabase.h
  SiloDatabaseFactory.h
  StartupShutdownManager.h
  Statistic.h
  Statistician.h
  StatTransaction.h
  Timer.h
  TimerManager.h
  Tracer.h
  Transaction.h
  Utilities.h)

set_source_files_properties(
  Array.C
  MathUtilities.C
  AsyncCommPeer.C
  PROPERTIES HEADER_FILE_ONLY TRUE)

set (tbox_sources
  ArraySpecial.C
  AsyncCommGroup.C
  AsyncCommStage.C
  BalancedDepthFirstTree.C
  BreadthFirstRankTree.C
  CenteredRankTree.C
  Clock.C
  CommGraphWriter.C
  Database.C
  DatabaseBox.C
  DatabaseFactory.C
  Dimension.C
  Grammar.C
  HDFDatabase.C
  HDFDatabaseFactory.C
  IEEE.C
  InputManager.C
  Logger.C
  MathUtilitiesSpecial.C
  MemoryDatabase.C
  MemoryDatabaseFactory.C
  MemoryUtilities.C
  MessageStream.C
  NullDatabase.C
  PIO.C
  ParallelBuffer.C
  Parser.C
  RankGroup.C
  RankTreeStrategy.C
  ReferenceCounter.C
  RestartManager.C
  SAMRAIManager.C
  SAMRAI_MPI.C
  Scanner.C
  Schedule.C
  Serializable.C
  SiloDatabase.C
  SiloDatabaseFactory.C
  StartupShutdownManager.C
  StatTransaction.C
  Statistic.C
  Statistician.C
  Timer.C
  TimerManager.C
  Tracer.C
  Transaction.C
  Utilities.C)

if (HAVE_HDF5)
    list(APPEND tbox_depends_on hdf5)
endif ()

if (ENABLE_MPI)
  list(APPEND tbox_depends_on mpi)
endif ()

blt_add_library( NAME SAMRAI_tbox
  SOURCES ${tbox_sources}
  HEADERS ${tbox_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/tbox
  DEPENDS_ON ${tbox_depends_on})

target_include_directories( SAMRAI_tbox
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SAMRAI_tbox
  EXPORT SAMRAITargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
