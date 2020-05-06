set(vector_depends ${SAMRAI_LIBRARIES})

if (ENABLE_OPENMP)
  set(vector_depends ${vector_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(vector_depends ${vector_depends} cuda) 
endif ()


blt_add_executable(
  NAME kvtest
  SOURCES kvtest.C
  DEPENDS_ON ${vector_depends})

blt_add_test(
  NAME kvtest
  COMMAND kvtest)

blt_add_executable(
  NAME pvtest
  SOURCES pvtest.C
  DEPENDS_ON ${vector_depends})

blt_add_test(
  NAME pvtest
  COMMAND pvtest)
