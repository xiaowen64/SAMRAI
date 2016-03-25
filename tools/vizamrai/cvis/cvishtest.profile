#!/bin/sh
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/cvis/cvishtest.profile $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
