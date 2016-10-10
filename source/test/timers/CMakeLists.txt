set ( timer_sources
        Foo.C
        Foo.h
        main_example.C
 
        main_timer.C)

set ( stats_sources
        main_stats.C

set ( example_sources
  main_example.C)

target_add_definitions (TESTING=1)

blt_add_executable(
  NAME timertest
)

blt_add_executable(
  NAME statstest
)

blt_add_executable(
  NAME example
  SOURCES ${exmaple_sources})
