c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 4d patch data
c               on a regular Cartesian mesh.
c
include(geom_m4cartcoarsenops4d.i)dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub4d(
cart_wgtavg_op_cell_4d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot4d(
cart_wgtavg_op_cell_4d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx4d(
cart_wgtavg_op_cell_4d(`double complex')dnl
