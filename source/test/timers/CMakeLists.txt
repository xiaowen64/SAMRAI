set ( timer_sources
        Foo.C
        main_timer.C)

set ( stats_sources
        main_stats.C)

set ( example_sources
  main_example.C)

blt_add_executable(
  NAME timer_timertest
  SOURCES ${timer_sources}
  DEPENDS_ON
    SAMRAI_tbox)

blt_add_executable(
  NAME timer_statstest
  SOURCES ${stats_sources}
  DEPENDS_ON
    SAMRAI_tbox)

blt_add_executable(
  NAME timer_example
  SOURCES ${example_sources}
  DEPENDS_ON
    SAMRAI_tbox)


target_compile_definitions(timer_timertest PUBLIC TESTING=1)
target_compile_definitions(timer_statstest PUBLIC TESTING=1)
target_compile_definitions(timer_example PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME timer_timertest
  EXECUTABLE timer_timertest
  INPUTS ${test_inputs}
  PARALLEL TRUE)

blt_add_test(
  NAME timer_statstest
  COMMAND timer_statstest)
