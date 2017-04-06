set (facadaptive_sources
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

set (facadaptive_depends_on
  testlib
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_mesh
  SAMRAI_solv
  SAMRAI_appu
  SAMRAI_xfer
  SAMRAI_geom
  SAMRAI_tbox)

blt_add_executable(
  NAME facadaptive
  SOURCES ${facadaptive_sources}
  DEPENDS_ON
    ${facadaptive_depends_on})

target_compile_definitions(facadaptive PUBLIC TESTING=1)
