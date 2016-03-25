#!/bin/sh
##
## File:        cvistest
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Script used to invoke the CASC Visualization toolkit
##              used for testing out of development directory
##

LOCAL=/home/ssmith/local/profile
LD_LIBRARY_PATH=$LOCAL/lib:$LD_LIBRARY_PATH
TCLLIBPATH=.
export LD_LIBRARY_PATH TCLLIBPATH
echo opengl with profiling
VTK_RENDERER=OpenGL
export VTK_RENDERER
$LOCAL/bin/vtk $*
