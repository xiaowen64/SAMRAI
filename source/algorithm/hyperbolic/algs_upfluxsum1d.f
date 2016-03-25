c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/hyperbolic/algs_upfluxsum1d.f $
c  Package:     SAMRAI algorithms
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: F77 routines for updating 1d flux sums from fluxes.
c
c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/hyperbolic/algs_upfluxsum1d.f $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 1d arrays in FORTRAN routines.
c
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
      subroutine upfluxsum1d(
     &  ilo0,ihi0,
     &  flxgc0,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ihi0,
     &  flxgc0,
     &  iface
      double precision
     &  flux(ilo0-flxgc0:ihi0+1+flxgc0),
     &  fluxsum(1)
      integer ie0
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie0 = ilo0
      else
        ie0 = ihi0+1
      endif 
 
      fluxsum(1)=fluxsum(1)+flux(ie0)
c
      return
      end
