set (mpi_on_sources
  main-runtime-on.C
  mpi-interface-tests.C)

set (mpi_off_sources
  main-runtime-off.C
  mpi-interface-tests.C)

set (mpi_disabled_sources
  main-runtime-disabled.C
  mpi-interface-tests.C)

blt_add_executable(
  NAME mpi_on
  SOURCES ${mpi_on_sources}
  DEPENDS_ON
    SAMRAI_tbox)

blt_add_executable(
  NAME mpi_off
  SOURCES ${mpi_off_sources}
  DEPENDS_ON
    SAMRAI_tbox)

blt_add_executable(
  NAME mpi_disabled
  SOURCES ${mpi_disabled_sources}
  DEPENDS_ON
    SAMRAI_tbox)

target_compile_definitions(mpi_on PUBLIC TESTING=1)
target_compile_definitions(mpi_off PUBLIC TESTING=1)
target_compile_definitions(mpi_disabled PUBLIC TESTING=1)
