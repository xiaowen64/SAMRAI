c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/fortran/pdat_dbugfort2d.m4 $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1704 $
c  Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
c  Description: F77 routines for debugging 2d patch data types.
c
include(pdat_dbugstuff2d.i)dnl
c
c***********************************************************************
c Debugging routines for 2d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub2d(
pdat_debug_cell_2d(`double precision')dnl
c
      subroutine dbugcellflot2d(
pdat_debug_cell_2d(`real')dnl
c
      subroutine dbugcellcplx2d(
pdat_debug_cell_2d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 2d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub2d0(
pdat_debug_face_2d(`double precision',0,1)dnl
c
      subroutine dbugfacedoub2d1(
pdat_debug_face_2d(`double precision',1,0)dnl
c
      subroutine dbugfaceflot2d0(
pdat_debug_face_2d(`real',0,1)dnl
c
      subroutine dbugfaceflot2d1(
pdat_debug_face_2d(`real',1,0)dnl
c
      subroutine dbugfacecplx2d0(
pdat_debug_face_2d(`double complex',0,1)dnl
c
      subroutine dbugfacecplx2d1(
pdat_debug_face_2d(`double complex',1,0)dnl
c
c
c***********************************************************************
c Debugging routines for 2d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub2d(
pdat_debug_node_2d(`double precision')dnl
c
      subroutine dbugnodeflot2d(
pdat_debug_node_2d(`real')dnl
c
      subroutine dbugnodecplx2d(
pdat_debug_node_2d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 2d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub2d0(
pdat_debug_outerface_2d(`double precision',0,1)dnl
c
      subroutine dbugoutfacedoub2d1(
pdat_debug_outerface_2d(`double precision',1,0)dnl
c
      subroutine dbugoutfaceflot2d0(
pdat_debug_outerface_2d(`real',0,1)dnl
c
      subroutine dbugoutfaceflot2d1(
pdat_debug_outerface_2d(`real',1,0)dnl
c
      subroutine dbugoutfacecplx2d0(
pdat_debug_outerface_2d(`double complex',0,1)dnl
c
      subroutine dbugoutfacecplx2d1(
pdat_debug_outerface_2d(`double complex',1,0)dnl
c
c
c***********************************************************************
c Debugging routines for 2d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub2d0(
pdat_debug_outerside_2d(`double precision',0,1)dnl
c
      subroutine dbugoutsidedoub2d1(
pdat_debug_outerside_2d(`double precision',1,0)dnl
c
      subroutine dbugoutsideflot2d0(
pdat_debug_outerside_2d(`real',0,1)dnl
c
      subroutine dbugoutsideflot2d1(
pdat_debug_outerside_2d(`real',1,0)dnl
c
      subroutine dbugoutsidecplx2d0(
pdat_debug_outerside_2d(`double complex',0,1)dnl
c
      subroutine dbugoutsidecplx2d1(
pdat_debug_outerside_2d(`double complex',1,0)dnl
c
c
c***********************************************************************
c Debugging routines for 2d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub2d0(
pdat_debug_side_2d(`double precision',0,1)dnl
c
      subroutine dbugsidedoub2d1(
pdat_debug_side_2d(`double precision',1,0)dnl
c
      subroutine dbugsideflot2d0(
pdat_debug_side_2d(`real',0,1)dnl
c
      subroutine dbugsideflot2d1(
pdat_debug_side_2d(`real',1,0)dnl
c
      subroutine dbugsidecplx2d0(
pdat_debug_side_2d(`double complex',0,1)dnl
c
      subroutine dbugsidecplx2d1(
pdat_debug_side_2d(`double complex',1,0)dnl
c
c
c***********************************************************************
c Debugging routines for 2d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub2d0(
pdat_debug_edge_2d(`double precision',0,1)dnl
c
      subroutine dbugedgedoub2d1(
pdat_debug_edge_2d(`double precision',1,0)dnl
c
      subroutine dbugedgeflot2d0(
pdat_debug_edge_2d(`real',0,1)dnl
c
      subroutine dbugedgeflot2d1(
pdat_debug_edge_2d(`real',1,0)dnl
c
      subroutine dbugedgecplx2d0(
pdat_debug_edge_2d(`double complex',0,1)dnl
c
      subroutine dbugedgecplx2d1(
pdat_debug_edge_2d(`double complex',1,0)dnl
c
