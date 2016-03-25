CMAKE_MINIMUM_REQUIRED(VERSION 1.4)
PROJECT(VIZAMRAI)

#-----------------------------------------------------------------------------
# Platform configuration tests.
INCLUDE(${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

SET (LIBRARY_OUTPUT_PATH ${VIZAMRAI_BINARY_DIR}/bin/ CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${VIZAMRAI_BINARY_DIR}/bin/ CACHE PATH "Single output directory for building all executables.")

SET (VIZAMRAI_CAN_BUILD 1)

# If using Windows need to link to an XDR library
IF (WIN32)
   SUBDIRS(vtkcvis ONCRPC Wrapping scripts)
   SUBDIR_DEPENDS(vtkcvis ONCRPC)
ELSE (WIN32)
   SUBDIRS(vtkcvis Wrapping scripts)
ENDIF (WIN32)

# include module to find vtk
INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
# if vtk was found include the use vtk file
IF (USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
ELSE(USE_VTK_FILE)
  SET (VIZAMRAI_CAN_BUILD 0)
  MESSAGE("vizamrai requires VTK to build. Please set VTK_BINARY_DIR to directory where vtkConfigure.h and UseVTK.cmake are found.")
ENDIF (USE_VTK_FILE)

IF (VIZAMRAI_CAN_BUILD)
 # Should add more checks here to make sure VTK is built correctly
 IF ( NOT VTK_WRAP_TCL )
   MESSAGE("vizamrai requires VTK to be built with Tcl bindings. Please rebuild VTK with Tcl bindings and re-run CMake")
   SET (VIZAMRAI_CAN_BUILD 0)
 ENDIF ( NOT VTK_WRAP_TCL )
ENDIF (VIZAMRAI_CAN_BUILD)

IF (VIZAMRAI_CAN_BUILD)
 OPTION(BUILD_SHARED_LIBS "Build vizamrai with shared libraries." ${VTK_BUILD_SHARED_LIBS})

 IF (TK_INCLUDE_PATH)
   INCLUDE_DIRECTORIES(${TK_INCLUDE_PATH})
 ENDIF (TK_INCLUDE_PATH)

 INCLUDE_DIRECTORIES(${VIZAMRAI_BINARY_DIR}/vtkcvis 
                     ${VIZAMRAI_SOURCE_DIR}/vtkcvis)

 IF (WIN32)
    INCLUDE_DIRECTORIES(${VIZAMRAI_SOURCE_DIR}/ONCRPC)
 ENDIF (WIN32)

 LINK_DIRECTORIES(${VIZAMRAI_BINARY_DIR}/vizcvis)

 IF (UNIX)
    LINK_LIBRARIES(${THREAD_LIBS} ${DL_LIBS} -lm)
 ENDIF (UNIX)

 FIND_FILE (VTK_WRAP_HINTS hints ${VIZAMRAI_SOURCE_DIR}/Wrapping)


#***************************************************************************
# This is modeled after FindTclsh.cmake script 
#
# This module finds if ITCLSH is installed and determines where the
# include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#
#  ITCL_ITCLSH        = the full path to the tclsh binary (tcl tcl80 etc)
#

# In cygwin, look for the cygwin version first.  Don't look for it later to
# avoid finding the cygwin version on a Win32 build.
IF(WIN32)
  IF(UNIX)
    FIND_PROGRAM(ITCL_ITCLSH cygitclsh80)
  ENDIF(UNIX)
ENDIF(WIN32)

FIND_PROGRAM(ITCL_ITCLSH
  NAMES itclsh
  itclsh84 itclsh8.4
  itclsh83 itclsh8.3
  itclsh82 itclsh8.2
  itclsh80 itclsh8.0
)

IF (WIN32)
  MARK_AS_ADVANCED(
    ITCL_ITCLSH
    )
ENDIF(WIN32)

#***************************************************************************
# This is modeled after FindTclsh.cmake script 
#
# This module finds if ITCLSH is installed and determines where the
# include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#
#  ITK_ITKWISH        = the full path to the tclsh binary (tcl tcl80 etc)
#

# In cygwin, look for the cygwin version first.  Don't look for it later to
# avoid finding the cygwin version on a Win32 build.
IF(WIN32)
  IF(UNIX)
    FIND_PROGRAM(ITK_ITKWISH cygitkwish80)
  ENDIF(UNIX)
ENDIF(WIN32)

FIND_PROGRAM(ITK_ITKWISH
  NAMES itkwish
  itkwish84 itkwish8.4
  itkwish83 itkwish8.3
  itkwish82 itkwish8.2
  itkwish80 itkwish8.0
)

IF (WIN32)
  MARK_AS_ADVANCED(
    ITCL_ITKWISH
    )
ENDIF(WIN32)
#****************************************************************************

# This is modeled after FindTCL.cmake
#
# This module finds if ITCL is installed and determines where the
# include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#
#  ITCL_LIBRARY            = full path to the Tcl library (tcl tcl80 etc)
#  ITCL_LIBRARY_DEBUG      = full path to the Tcl library (debug)
#  ITCL_STUB_LIBRARY       = full path to the Tcl stub library
#  ITCL_STUB_LIBRARY_DEBUG = full path to the Tcl stub library (debug)
#  ITCL_INCLUDE_PATH       = path to where tcl.h can be found
#  ITCL_ITCLSH              = full path to the tclsh binary (tcl tcl80 etc)
#  ITK_LIBRARY             = full path to the Tk library (tk tk80 etc)
#  ITK_LIBRARY_DEBUG       = full path to the Tk library (debug)
#  ITK_STUB_LIBRARY        = full path to the Tk stub library
#  ITK_STUB_LIBRARY_DEBUG  = full path to the Tk stub library (debug)
#  ITK_INCLUDE_PATH        = path to where tk.h can be found
#  ITK_INTERNAL_PATH       = path to where tkWinInt.h can be found
#  ITK_ITKWISH                = full path to the itkwish binary (itkwish itkwish80 etc)
#

GET_FILENAME_COMPONENT(ITCL_ITCLSH_PATH ${ITCL_ITCLSH} PATH)

GET_FILENAME_COMPONENT(ITK_ITKWISH_PATH ${ITK_ITKWISH} PATH)

SET (ITCLITK_POSSIBLE_LIB_PATHS
  "${ITCL_ITCLSH_PATH}/../lib"
  "${ITK_ITKWISH_PATH}/../lib"
  /usr/lib 
  /usr/local/lib
  "C:/Program Files/Tcl/lib" 
  "C:/Tcl/lib" 
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.4;Root]/lib
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.3;Root]/lib
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.2;Root]/lib
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.0;Root]/lib
)

FIND_LIBRARY(ITCL_LIBRARY
  NAMES itcl itcl84 itcl8.4 itcl83 itcl8.3 itcl82 itcl8.2 itcl80 itcl8.0
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITCL_LIBRARY_DEBUG
  NAMES itcld itcl84d itcl8.4d itcl83d itcl8.3d itcl82d itcl8.2d itcl80d itcl8.0d
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITCL_STUB_LIBRARY
  NAMES itclstub itclstub84 itclstub8.4 itclstub83 itclstub8.3 itclstub82 itclstub8.2 itclstub80 itclstub8.0
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITCL_STUB_LIBRARY_DEBUG
  NAMES itclstubd itclstub84d itclstub8.4d itclstub83d itclstub8.3d itclstub82d itclstub8.2d itclstub80d itclstub8.0d
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITK_LIBRARY 
  NAMES itk itk84 itk8.4 itk83 itk8.3 itk82 itk8.2 itk80 itk8.0
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITK_LIBRARY_DEBUG
  NAMES itkd itk84d itk8.4d itk83d itk8.3d itk82d itk8.2d itk80d itk8.0d
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITK_STUB_LIBRARY 
  NAMES itkstub itkstub84 itkstub8.4 itkstub83 itkstub8.3 itkstub82 itkstub8.2 itkstub80 itkstub8.0
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

FIND_LIBRARY(ITK_STUB_LIBRARY_DEBUG
  NAMES itkstubd itkstub84d itkstub8.4d itkstub83d itkstub8.3d itkstub82d itkstub8.2d itkstub80d itkstub8.0d
  PATHS ${ITCLITK_POSSIBLE_LIB_PATHS}
)

SET (ITCLITK_POSSIBLE_INCLUDE_PATHS
  "${ITCL_ITCLSH_PATH}/../include"
  "${ITK_ITKWISH_PATH}/../include"
  /usr/include 
  /usr/local/include
  /usr/include/tcl8.4
  /usr/include/tcl8.3
  /usr/include/tcl8.2
  /usr/include/tcl8.0
  "C:/Program Files/Tcl/include" 
  "C:/Tcl/include" 
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.4;Root]/include
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.3;Root]/include
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.2;Root]/include
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Scriptics\\Tcl\\8.0;Root]/include
)

FIND_PATH(ITCL_INCLUDE_PATH tcl.h 
  ${ITCLITK_POSSIBLE_INCLUDE_PATHS}
)

FIND_PATH(ITK_INCLUDE_PATH tk.h 
  ${ITCLITK_POSSIBLE_INCLUDE_PATHS}
)

IF (WIN32)
  FIND_PATH(ITK_INTERNAL_PATH tkWinInt.h
    ${ITCLITK_POSSIBLE_INCLUDE_PATHS}
  )
  MARK_AS_ADVANCED(ITK_INTERNAL_PATH)

  MARK_AS_ADVANCED(
    ITCL_ITCLSH_PATH
    ITK_ITKWISH_PATH
    ITCL_INCLUDE_PATH
    ITK_INCLUDE_PATH
    ITCL_LIBRARY
    ITCL_LIBRARY_DEBUG
    ITCL_STUB_LIBRARY
    ITCL_STUB_LIBRARY_DEBUG
    ITK_LIBRARY  
    ITK_LIBRARY_DEBUG
    ITK_STUB_LIBRARY  
    ITK_STUB_LIBRARY_DEBUG
    )
ENDIF(WIN32)

#***************************************************************************


ENDIF (VIZAMRAI_CAN_BUILD)




