define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine error(
     & lo0, hi0, 
     & u, w,
     & lambda,
     & xlo, xhi, dx,
     & t,
     & maxerror,
     & l2error )

c  Compute error in solution at time t.

      implicit none

      integer lo0

      integer hi0

      double precision lambda
      double precision t
      double precision l2error
      double precision maxerror

      double precision u(CELL1d(lo,hi,0))
      double precision w(CELL1d(lo,hi,0))

      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)
      double precision dx(0:NDIM-1)

      integer i

      double precision diff
      double precision error_i
      double precision localerror

      double precision xi, xterm

      double precision zero,       half,       one
      parameter      ( zero=0.0d0, half=0.5d0, one=1.0d0 )

      intrinsic abs, max
      double precision abs, max

      localerror = zero
      l2error = zero
      xi = xlo(0) + dx(0)*half
      do i = lo0, hi0
         xterm = xi*(one - xi)
         error_i = abs(u(i) - t*xterm)
         l2error = l2error + w(i)*error_i**2
         diff = w(i)*error_i
         localerror = max(localerror, diff)
         xi = xi + dx(0)
      end do

c  Since the weight was used to mask coarse cells that have been refined,
c  the error in those that aren't covered are scaled by the cell volume.

      localerror = localerror/dx(0)
      l2error = sqrt(l2error)

      maxerror = max(localerror, maxerror)

      return
      end
