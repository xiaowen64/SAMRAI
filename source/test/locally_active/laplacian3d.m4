define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

c
c File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/locally_active/laplacian3d.m4 $
c Package:     SAMRAI application
c Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
c Revision:    $LastChangedRevision: 1704 $
c Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
c Description: Fortran code for locally active tests.
c

c***********************************************************************
c***********************************************************************
      subroutine laplacian(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0, ngc1, ngc2,
     &  dx,
     &  func, soln)
c***********************************************************************
      implicit none
c***********************************************************************
c inputs:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngc0,ngc1,ngc2
      double precision dx(0:2)
      double precision 
     &  func(CELL3dVECG(ifirst,ilast,ngc)),
     &  soln(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision d2ux, d2uy, d2uz

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
 
               d2ux = func(ic0+1,ic1,ic2) - 2.0*func(ic0,ic1,ic2) +        
     &                func(ic0-1,ic1,ic2)
               d2uy = func(ic0,ic1+1,ic2) - 2.0*func(ic0,ic1,ic2) + 
     &                func(ic0,ic1-1,ic2)
               d2uz = func(ic0,ic1,ic2+1) - 2.0*func(ic0,ic1,ic2) + 
     &                func(ic0,ic1,ic2-1)
               soln(ic0,ic1,ic2) = d2ux/dx(0)**2 + d2uy/dx(1)**2 +
     &                d2uz/dx(2)**2

            enddo
         enddo
      enddo

      return
      end



c***********************************************************************
c***********************************************************************
      subroutine init(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  xlo,dx, 
     &  centroid, alpha,
     &  func)
c***********************************************************************
      implicit none
      double precision half, smallr
      parameter(half=0.5d0)
      parameter(smallr=1.0e-6)
c***********************************************************************
c inputs:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngc0, ngc1, ngc2
      double precision xlo(0:NDIM-1), dx(0:NDIM-1)
      double precision centroid(0:NDIM-1), alpha
      double precision 
     &  func(CELL3dVECG(ifirst,ilast,ngc))
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision xc(0:NDIM-1), radius, x0, x1, x2

c
c  Initialize patch interiors.
c
      do ic2=ifirst2,ilast2
         xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
         x2 = xc(2)-centroid(2)

         do ic1=ifirst1,ilast1
            xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
            x1 = xc(1)-centroid(1)

            do ic0=ifirst0,ilast0
               xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
               x0 = xc(0)-centroid(0)

               radius = sqrt(max(smallr,(x0**2 + x1**2 + x2**2)))
               func(ic0,ic1,ic2) = exp(-radius*alpha)
        
            end do
         end do
      end do

      return
      end


c***********************************************************************
c***********************************************************************
      subroutine bdry(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  bfirst0,blast0,bfirst1,blast1,bfirst2,blast2,
     &  xlo,dx, 
     &  centroid, alpha,
     &  func)
c***********************************************************************
      implicit none
      double precision half, smallr
      parameter(half=0.5d0)
      parameter(smallr=1.0e-6)
c***********************************************************************
c inputs:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngc0,ngc1,ngc2
      integer bfirst0,blast0,bfirst1,blast1,bfirst2,blast2
      double precision xlo(0:NDIM-1), dx(0:NDIM-1)
      double precision centroid(0:NDIM-1), alpha
      double precision 
     &  func(CELL3dVECG(ifirst,ilast,ngc))
c***********************************************************************

      integer ic0,ic1,ic2
      double precision xc(0:NDIM-1), radius, x0, x1, x2

c
c  Initialize  fill box [(bfirst),(blast)] values.
c
      do ic2=bfirst2,blast2
         xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
         x2 = xc(2)-centroid(2)

         do ic1=bfirst1,blast1
            xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
            x1 = xc(1)-centroid(1)

            do ic0=bfirst0,blast0
               xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
               x0 = xc(0)-centroid(0)

               radius = dsqrt(max(smallr,(x0**2 + x1**2 + x2**2)))
               func(ic0,ic1,ic2) = exp(-radius*alpha)

            end do
        
         end do

      end do

      return
      end


