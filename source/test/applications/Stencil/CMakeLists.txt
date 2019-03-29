set (stencil_sources
  Stencil.C
  main.C)

set(stencil_depends ${SAMRAI_LIBRARIES})

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(stencil_depends ${stencil_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(stencil_depends ${stencil_depends} cuda)
endif ()

blt_add_executable(
  NAME stencil
  SOURCES ${stencil_sources}
  DEPENDS_ON ${stencil_depends})

target_include_directories( stencil
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/applications/Stencil)

set_target_properties(stencil PROPERTIES LINKER_LANGUAGE CXX)
