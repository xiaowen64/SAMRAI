c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/cartesian/operators/fortran/geom_cartcoarsen1d.m4 $
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: FORTRAN routines for spatial coarsening of 1d patch data
c               on a regular Cartesian mesh.
c
include(geom_m4cartcoarsenops1d.i)dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub1d(
cart_wgtavg_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot1d(
cart_wgtavg_op_cell_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx1d(
cart_wgtavg_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub1d(
cart_wgtavg_op_edge_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot1d(
cart_wgtavg_op_edge_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx1d(
cart_wgtavg_op_edge_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub1d(
cart_wgtavg_op_face_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot1d(
cart_wgtavg_op_face_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx1d(
cart_wgtavg_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub1d(
cart_wgtavg_op_outerface_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot1d(
cart_wgtavg_op_outerface_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx1d(
cart_wgtavg_op_outerface_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub1d(
cart_wgtavg_op_outerside_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub1d(
cart_wgtavg_op_side_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot1d(
cart_wgtavg_op_side_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx1d(
cart_wgtavg_op_side_1d(`double complex')dnl
c
