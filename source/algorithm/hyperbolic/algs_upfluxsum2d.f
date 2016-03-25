c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/hyperbolic/algs_upfluxsum2d.f $
c  Package:     SAMRAI algorithms
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for updating 2d flux sums from fluxes.
c
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/hyperbolic/algs_upfluxsum2d.f $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
      subroutine upfluxsumface2d0(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1),
     &  fluxsum(ilo1:ihi1)
      integer ie0,ic1
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie0 = ilo0
      else 
        ie0 = ihi0+1
      endif 
 
      do ic1=ilo1,ihi1
         fluxsum(ic1)=fluxsum(ic1)+flux(ie0,ic1)
      enddo
c
      return
      end
c
      subroutine upfluxsumface2d1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      double precision
     &  flux(ilo1-flxgc1:ihi1+1+flxgc1,
     &          ilo0-flxgc0:ihi0+flxgc0),
     &  fluxsum(ilo0:ihi0)
      integer ie1,ic0
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie1 = ilo1
      else 
        ie1 = ihi1+1
      endif 
 
      do ic0=ilo0,ihi0
         fluxsum(ic0)=fluxsum(ic0)+flux(ie1,ic0)
      enddo
c
      return
      end
c
      subroutine upfluxsumside2d0(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0,
     &          ilo1-flxgc1:ihi1+flxgc1),
     &  fluxsum(ilo1:ihi1)
      integer ic0,ic1
c
c***********************************************************************
c
      if (iface.eq.0) then
        ic0 = ilo0
      else 
        ic0 = ihi0+1
      endif 
 
      do ic1=ilo1,ihi1
         fluxsum(ic1)=fluxsum(ic1)+flux(ic0,ic1)
      enddo
c
      return
      end
c
      subroutine upfluxsumside2d1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+flxgc0,
     &          ilo1-flxgc1:ihi1+1+flxgc1),
     &  fluxsum(ilo0:ihi0)
      integer ic1,ic0
c
c***********************************************************************
c
      if (iface.eq.0) then
        ic1 = ilo1
      else 
        ic1 = ihi1+1
      endif 
 
      do ic0=ilo0,ihi0
         fluxsum(ic0)=fluxsum(ic0)+flux(ic0,ic1)
      enddo
c
      return
      end
c
