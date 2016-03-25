c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/fortran/pdat_dbugstuff1d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 1d patchdata debugging routines
c
define(SAMRAI_FORTDIR,.)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
define(pdat_debug_subroutine_head_1d,`dnl
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
')dnl
c
define(pdat_debug_body_1d,`dnl

      do $1=$2
         write(6,*) "array[",$1,"] = ",array($1)
      enddo

      call flush(6)
      return
      end
')dnl
define(pdat_debug_cell_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(CELL1d(fi,la,ng))
      integer ic0
c     =============================================================
pdat_debug_body_1d(`ic0',`ibeg0,iend0')dnl
')dnl
define(pdat_debug_face_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(FACE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_node_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(NODE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_outerface_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(OUTERFACE1d(fi,la,ng))
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
')dnl
define(pdat_debug_outerside_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(OUTERSIDE1d(fi,la,ng))
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
')dnl
define(pdat_debug_side_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(SIDE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_edge_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(EDGE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0')dnl
')dnl
