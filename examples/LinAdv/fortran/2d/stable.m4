define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine stabledt(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  advecspeed,uval,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ngc0,ngc1
c
      REAL  
     &  advecspeed(0:NDIM-1),
     &  uval(CELL2dVECG(ifirst,ilast,ngc))
c    
      REAL maxspeed(0:NDIM-1)
c
      maxspeed(0)=zero
      maxspeed(1)=zero

      maxspeed(0) = max(maxspeed(0), abs(advecspeed(0)))
      maxspeed(1) = max(maxspeed(1), abs(advecspeed(1)))
      stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
      return
      end 
