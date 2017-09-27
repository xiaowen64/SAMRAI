set (geom_headers
        CartesianCellComplexConservativeLinearRefine.h
        CartesianCellComplexLinearRefine.h
        CartesianCellComplexWeightedAverage.h
        CartesianCellDoubleConservativeLinearRefine.h
        CartesianCellDoubleLinearRefine.h
        CartesianCellDoubleWeightedAverage.h
        CartesianCellFloatConservativeLinearRefine.h
        CartesianCellFloatLinearRefine.h
        CartesianCellFloatWeightedAverage.h
        CartesianEdgeComplexWeightedAverage.h
        CartesianEdgeDoubleConservativeLinearRefine.h
        CartesianEdgeDoubleWeightedAverage.h
        CartesianEdgeFloatConservativeLinearRefine.h
        CartesianEdgeFloatWeightedAverage.h
        CartesianFaceComplexWeightedAverage.h
        CartesianFaceDoubleConservativeLinearRefine.h
        CartesianFaceDoubleWeightedAverage.h
        CartesianFaceFloatConservativeLinearRefine.h
        CartesianFaceFloatWeightedAverage.h
        CartesianGridGeometry.h
        CartesianNodeComplexLinearRefine.h
        CartesianNodeDoubleLinearRefine.h
        CartesianNodeFloatLinearRefine.h
        CartesianOuterfaceComplexWeightedAverage.h
        CartesianOuterfaceDoubleWeightedAverage.h
        CartesianOuterfaceFloatWeightedAverage.h
        CartesianOutersideDoubleWeightedAverage.h
        CartesianPatchGeometry.h
        CartesianSideComplexWeightedAverage.h
        CartesianSideDoubleConservativeLinearRefine.h
        CartesianSideDoubleWeightedAverage.h
        CartesianSideFloatConservativeLinearRefine.h
        CartesianSideFloatWeightedAverage.h
        GridGeometry.h)

set (geom_sources
  CartesianCellComplexConservativeLinearRefine.C
  CartesianCellComplexLinearRefine.C
  CartesianCellComplexWeightedAverage.C
  CartesianCellDoubleConservativeLinearRefine.C
  CartesianCellDoubleLinearRefine.C
  CartesianCellDoubleWeightedAverage.C
  CartesianCellFloatConservativeLinearRefine.C
  CartesianCellFloatLinearRefine.C
  CartesianCellFloatWeightedAverage.C
  CartesianEdgeComplexWeightedAverage.C
  CartesianEdgeDoubleConservativeLinearRefine.C
  CartesianEdgeDoubleWeightedAverage.C
  CartesianEdgeFloatConservativeLinearRefine.C
  CartesianEdgeFloatWeightedAverage.C
  CartesianFaceComplexWeightedAverage.C
  CartesianFaceDoubleConservativeLinearRefine.C
  CartesianFaceDoubleWeightedAverage.C
  CartesianFaceFloatConservativeLinearRefine.C
  CartesianFaceFloatWeightedAverage.C
  CartesianGridGeometry.C
  CartesianNodeComplexLinearRefine.C
  CartesianNodeDoubleLinearRefine.C
  CartesianNodeFloatLinearRefine.C
  CartesianOuterfaceComplexWeightedAverage.C
  CartesianOuterfaceDoubleWeightedAverage.C
  CartesianOuterfaceFloatWeightedAverage.C
  CartesianOutersideDoubleWeightedAverage.C
  CartesianPatchGeometry.C
  CartesianSideComplexWeightedAverage.C
  CartesianSideDoubleConservativeLinearRefine.C
  CartesianSideDoubleWeightedAverage.C
  CartesianSideFloatConservativeLinearRefine.C
  CartesianSideFloatWeightedAverage.C
  GridGeometry.C
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartcoarsen1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartcoarsen2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartcoarsen3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartcoarsen4d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartrefine1d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartrefine2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/geom_cartrefine3d.f)

set (geom_depends_on
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_tbox
  boost)

if (ENABLE_CUDA)
  set (geom_depends_on
    ${geom_depends_on}
    cuda)
endif ()

process_m4(NAME fortran/geom_cartcoarsen1d)
process_m4(NAME fortran/geom_cartcoarsen2d)
process_m4(NAME fortran/geom_cartcoarsen3d)
process_m4(NAME fortran/geom_cartcoarsen4d)
process_m4(NAME fortran/geom_cartrefine1d)
process_m4(NAME fortran/geom_cartrefine2d)
process_m4(NAME fortran/geom_cartrefine3d)

blt_add_library(
  NAME SAMRAI_geom
  SOURCES ${geom_sources}
  HEADERS ${geom_headers}
  HEADERS_OUTPUT_SUBDIR SAMRAI/geom
  DEPENDS_ON ${geom_depends_on})

install(TARGETS
  SAMRAI_geom
  RUNTIME DESTINATION lib
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

set_target_properties(SAMRAI_geom PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories( SAMRAI_geom
  PUBLIC ${PROJECT_SOURCE_DIR}/source)
