set (defaultlogger_sources
        defaultloggertest.C)

blt_add_executable(
        NAME defaultlogger
        SOURCES ${defaultlogger_sources}
        DEPENDS_ON
          SAMRAI_tbox)

set (userlogger_sources
        userloggertest.C)

blt_add_executable(
        NAME userlogger
        SOURCES ${userlogger_sources}
        DEPENDS_ON
        SAMRAI_tbox)