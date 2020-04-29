set (facadaptive_sources
  AdaptivePoisson.C
  GaussianFcn.C
  PoissonGaussianDiffcoefSolution.C
  PoissonGaussianSolution.C
  PoissonMultigaussianSolution.C
  PoissonPolynomialSolution.C
  PoissonSineSolution.C
  QuarticFcn.C
  SinusoidFcn.C
  main.C
  patchFcns.C
  setArrayData.C)

set(facadaptive_depends testlib)
set(facadaptive_depends ${facadaptive_depends} ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(facadaptive_depends ${facadaptive_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(facadaptive_depends ${facadaptive_depends} cuda)
endif ()

blt_add_executable(
  NAME facadaptive
  SOURCES ${facadaptive_sources}
  DEPENDS_ON
    ${facadaptive_depends})

target_compile_definitions(facadaptive PUBLIC TESTING=1)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)

samrai_add_tests(
  NAME facadaptive
  EXECUTABLE facadaptive
  INPUTS ${test_inputs}
  PARALLEL TRUE)
