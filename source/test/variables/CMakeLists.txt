
set(variables_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(variables_depends ${variables_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(variables_depends ${variables_depends} cuda)
endif ()

blt_add_executable(NAME variables
  SOURCES main-var_db.C
  DEPENDS_ON ${variables_depends})

target_compile_definitions(variables PUBLIC TESTING=1)

blt_add_test(
  NAME variables
  COMMAND variables)
