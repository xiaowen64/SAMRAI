# create the vtk executable
CONFIGURE_FILE(
  ${VIZAMRAI_SOURCE_DIR}/Wrapping/Tcl/tkAppInit.cxx.in
  ${VIZAMRAI_BINARY_DIR}/Wrapping/Tcl/tkAppInit.cxx)

SUBDIRS(Tcl)

IF (VTK_WRAP_TCL)
  INCLUDE (${CMAKE_ROOT}/Modules/FindTCL.cmake)
ENDIF (VTK_WRAP_TCL)

IF (TCL_LIBRARY)
  SET(LIBS ${LIBS} ${TCL_LIBRARY})
ELSE (TCL_LIBRARY)
  MESSAGE("VTK_WRAP_TCL is ON, but the Tcl library was not found.  Please set TCL_LIBRARY." "Error")
ENDIF (TCL_LIBRARY)


IF (VTK_USE_RENDERING)
  # add in the Tk values if found
  IF (TK_INCLUDE_PATH)
    INCLUDE_DIRECTORIES(${TK_INCLUDE_PATH})
    IF (APPLE)
      INCLUDE_DIRECTORIES(${VTK_SOURCE_DIR}/Rendering/tkOSXInternals)
    ENDIF (APPLE)
  ENDIF (TK_INCLUDE_PATH)
  IF (TK_XLIB_PATH)
    INCLUDE_DIRECTORIES(${TK_XLIB_PATH})
  ENDIF (TK_XLIB_PATH)

  SET(LIBS ${LIBS} vtkRenderingTCL vtkRendering)
ENDIF (VTK_USE_RENDERING)

IF (VTK_USE_PATENTED)
  SET(LIBS ${LIBS} vtkPatentedTCL vtkPatented)
ENDIF (VTK_USE_PATENTED)

IF (VTK_USE_HYBRID)
  SET(LIBS ${LIBS} vtkHybridTCL vtkHybrid)
ENDIF (VTK_USE_HYBRID)

IF (VTK_USE_PARALLEL)
  SET(LIBS ${LIBS} vtkParallelTCL vtkParallel)
ENDIF (VTK_USE_PARALLEL)

#
# link in the PNG library if it was found
#
ADD_EXECUTABLE(cvish ${VIZAMRAI_BINARY_DIR}/Wrapping/Tcl/tkAppInit.cxx)

IF (ITCL_LIBRARY)
SET(LIBS vtkcvisTCL vtkcvis ${ITCL_LIBRARY} ${LIBS})
ELSE (ITCL_LIBRARY)
  MESSAGE("The Incr Tcl library was not found.  Please set ITCL_LIBRARY." "Error")
ENDIF (ITCL_LIBRARY)

IF (ITK_LIBRARY)
SET(LIBS vtkcvisTCL vtkcvis ${ITK_LIBRARY} ${LIBS})
ELSE (ITK_LIBRARY)
  MESSAGE("The Incr Tcl library was not found.  Please set ITK_LIBRARY." "Error")
ENDIF (ITK_LIBRARY)

TARGET_LINK_LIBRARIES (cvish ${LIBS})

INSTALL_TARGETS(/bin cvish)