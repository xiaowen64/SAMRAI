
SET (cvis_SRCS
vtkSamraiStructuredPointsReader.cxx
vtkSamraiStructuredPointsVectorReader.cxx
vtkcvisScalarBarActor.cxx
vtkcvisLookupTable.cxx
vtkcvisHedgeHog.cxx
)

ADD_LIBRARY(vtkcvis ${cvis_SRCS})

IF (WIN32)
  LINK_LIBRARIES(oncrpc wsock32)
ENDIF (WIN32)

LINK_LIBRARIES(vtkRendering vtkHybrid vtkGraphics vtkFiltering vtkCommon )

IF (VTK_WRAP_TCL)
  VTK_WRAP_TCL(vtkcvisTCL cvisTCL_SRCS ${cvis_SRCS})
  ADD_LIBRARY(vtkcvisTCL ${cvisTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkcvisTCL vtkRenderingTCL vtkHybridTCL vtkGraphicsTCL vtkFilteringTCL  vtkCommonTCL  ${TCL_LIBRARY})
  INSTALL_TARGETS(/lib/vtk vtkcvisTCL)
ENDIF (VTK_WRAP_TCL)

INSTALL_TARGETS(/lib/vtk vtkcvis)
INSTALL_FILES(/include/vtk .h ${cvis_SRCS})

LINK_LIBRARIES(
vtkcvis
)

CONFIGURE_FILE(
  ${VIZAMRAI_SOURCE_DIR}/vtkcvis/vtkcvisConfigure.h.in 
  ${VIZAMRAI_BINARY_DIR}/vtkcvis/vtkcvisConfigure.h
)



