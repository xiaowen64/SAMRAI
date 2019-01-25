set (indx_dataops_sources
  indx_dataops.C
  SampleIndexData.C)

blt_add_executable(
  NAME indx_dataops
  SOURCES ${indx_dataops_sources}
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME cell_patchtest
  SOURCES cell_patchtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME cell_hiertest
  SOURCES cell_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME cell_hiertest3d
  SOURCES cell_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME edge_hiertest
  SOURCES edge_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME face_hiertest
  SOURCES face_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME side_hiertest
  SOURCES side_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME node_hiertest
  SOURCES node_hiertest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME cell_cplxtest
  SOURCES cell_cplxtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME edge_cplxtest
  SOURCES edge_cplxtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME node_cplxtest
  SOURCES node_cplxtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME face_cplxtest
  SOURCES face_cplxtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

blt_add_executable(
  NAME side_cplxtest
  SOURCES side_cplxtest.C
  DEPENDS_ON
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_pdat
    SAMRAI_math
    SAMRAI_tbox)

target_compile_definitions(indx_dataops PUBLIC DISPLAY TESTING=1)
target_compile_definitions(cell_patchtest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(cell_hiertest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(cell_hiertest3d PUBLIC DISPLAY TESTING=1)
target_compile_definitions(edge_hiertest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(face_hiertest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(side_hiertest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(node_hiertest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(cell_cplxtest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(edge_cplxtest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(node_cplxtest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(face_cplxtest PUBLIC DISPLAY TESTING=1)
target_compile_definitions(side_cplxtest PUBLIC DISPLAY TESTING=1)

target_include_directories(indx_dataops PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(cell_patchtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(cell_hiertest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(cell_hiertest3d PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(edge_hiertest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(face_hiertest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(side_hiertest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(node_hiertest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(cell_cplxtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(edge_cplxtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(node_cplxtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(face_cplxtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)
target_include_directories(side_cplxtest PUBLIC ${PROJECT_SOURCE_DIR}/source/test/dataops)

blt_add_test(
  NAME indx_dataops
  COMMAND indx_dataops)

blt_add_test(
  NAME cell_patchtest
  COMMAND cell_patchtest)

blt_add_test(
  NAME cell_hiertest
  COMMAND cell_hiertest)

blt_add_test(
  NAME cell_hiertest3d
  COMMAND cell_hiertest3d)

blt_add_test(
  NAME edge_hiertest
  COMMAND edge_hiertest)

blt_add_test(
  NAME face_hiertest
  COMMAND face_hiertest)

blt_add_test(
  NAME side_hiertest
  COMMAND side_hiertest)

blt_add_test(
  NAME node_hiertest
  COMMAND node_hiertest)

blt_add_test(
  NAME cell_cplxtest
  COMMAND cell_cplxtest)

blt_add_test(
  NAME edge_cplxtest
  COMMAND edge_cplxtest)

blt_add_test(
  NAME node_cplxtest
  COMMAND node_cplxtest)

blt_add_test(
  NAME face_cplxtest
  COMMAND face_cplxtest)

blt_add_test(
  NAME side_cplxtest
  COMMAND side_cplxtest)
