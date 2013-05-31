c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 4d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 4d Cartesian coarsen operators
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2013 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 4d arrays in FORTRAN routines.
c
c
c
c
c
c***********************************************************************
c Weighted averaging for 4d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub4d(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision dVf,dVc
      integer ic0,ic1,ic2,ic3,if0,if1,if2,if3,ir0,ir1,ir2,ir3
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)*dxf(3)
      dVc = dxc(0)*dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=zero
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ir0=0,ratio(0)-1
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ic1=ifirstc1,ilastc1
                           if1=ic1*ratio(1)+ir1
                           do ic0=ifirstc0,ilastc0
                              if0=ic0*ratio(0)+ir0
                              arrayc(ic0,ic1,ic2,ic3)=
     &                           arrayc(ic0,ic1,ic2,ic3)+
     &                           arrayf(if0,if1,if2,if3)*dVf
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=arrayc(ic0,ic1,ic2,ic3)/dVc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 4d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot4d(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision dVf,dVc
      integer ic0,ic1,ic2,ic3,if0,if1,if2,if3,ir0,ir1,ir2,ir3
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)*dxf(3)
      dVc = dxc(0)*dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=zero
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ir0=0,ratio(0)-1
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ic1=ifirstc1,ilastc1
                           if1=ic1*ratio(1)+ir1
                           do ic0=ifirstc0,ilastc0
                              if0=ic0*ratio(0)+ir0
                              arrayc(ic0,ic1,ic2,ic3)=
     &                           arrayc(ic0,ic1,ic2,ic3)+
     &                           arrayf(if0,if1,if2,if3)*dVf
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=arrayc(ic0,ic1,ic2,ic3)/dVc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 4d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx4d(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision dVf,dVc
      integer ic0,ic1,ic2,ic3,if0,if1,if2,if3,ir0,ir1,ir2,ir3
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)*dxf(3)
      dVc = dxc(0)*dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=cmplx(zero,zero)
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ir0=0,ratio(0)-1
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ic1=ifirstc1,ilastc1
                           if1=ic1*ratio(1)+ir1
                           do ic0=ifirstc0,ilastc0
                              if0=ic0*ratio(0)+ir0
                              arrayc(ic0,ic1,ic2,ic3)=
     &                           arrayc(ic0,ic1,ic2,ic3)+
     &                           arrayf(if0,if1,if2,if3)*dVf
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  arrayc(ic0,ic1,ic2,ic3)=arrayc(ic0,ic1,ic2,ic3)/dVc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 4d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub4d0(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision volf,volc
      integer ie0,ic1,ic2,ic3,if0,if1,if2,if3,ir1,ir2,ir3
c
c***********************************************************************
c
      volf=dxf(1)*dxf(2)*dxf(3)
      volc=dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=zero
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ic3=ifirstc3,ilastc3
                  if3=ic3*ratio(3)+ir3
                  do ic2=ifirstc2,ilastc2
                     if2=ic2*ratio(2)+ir2
                     do ic1=ifirstc1,ilastc1
                        if1=ic1*ratio(1)+ir1
                        do ie0=ifirstc0,ilastc0+1
                           if0=ie0*ratio(0)
                           arrayc(ie0,ic1,ic2,ic3)=
     &                        arrayc(ie0,ic1,ic2,ic3)
     &                           +arrayf(if0,if1,if2,if3)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=
     &               arrayc(ie0,ic1,ic2,ic3)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub4d1(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double precision
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo3:fihi3,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo3:cihi3,
     &          cilo0:cihi0)
      double precision volf,volc
      integer ie1,ic2,ic3,ic0,if1,if2,if3,if0,ir2,ir3,ir0
c
c***********************************************************************
c
      volf=dxf(2)*dxf(3)*dxf(0)
      volc=dxc(2)*dxc(3)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=zero
               enddo
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir3=0,ratio(3)-1
            do ir2=0,ratio(2)-1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ie1=ifirstc1,ilastc1+1
                           if1=ie1*ratio(1)
                           arrayc(ie1,ic2,ic3,ic0)=
     &                        arrayc(ie1,ic2,ic3,ic0)
     &                           +arrayf(if1,if2,if3,if0)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=
     &               arrayc(ie1,ic2,ic3,ic0)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub4d2(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double precision
     &  arrayf(filo2:fihi2+1,
     &          filo3:fihi3,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo3:cihi3,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision volf,volc
      integer ie2,ic3,ic0,ic1,if2,if3,if0,if1,ir3,ir0,ir1
c
c***********************************************************************
c
      volf=dxf(3)*dxf(0)*dxf(1)
      volc=dxc(3)*dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=zero
               enddo
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ir3=0,ratio(3)-1
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ic0=ifirstc0,ilastc0
                     if0=ic0*ratio(0)+ir0
                     do ic3=ifirstc3,ilastc3
                        if3=ic3*ratio(3)+ir3
                        do ie2=ifirstc2,ilastc2+1
                           if2=ie2*ratio(2)
                           arrayc(ie2,ic3,ic0,ic1)=
     &                        arrayc(ie2,ic3,ic0,ic1)
     &                           +arrayf(if2,if3,if0,if1)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=
     &               arrayc(ie2,ic3,ic0,ic1)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub4d3(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double precision
     &  arrayf(filo3:fihi3+1,
     &          filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo3:cihi3+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision volf,volc
      integer ie3,ic0,ic1,ic2,if3,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      volf=dxf(0)*dxf(1)*dxf(2)
      volc=dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=zero
               enddo
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        do ie3=ifirstc3,ilastc3+1
                           if3=ie3*ratio(3)
                           arrayc(ie3,ic0,ic1,ic2)=
     &                        arrayc(ie3,ic0,ic1,ic2)
     &                           +arrayf(if3,if0,if1,if2)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=
     &               arrayc(ie3,ic0,ic1,ic2)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 4d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot4d0(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision volf,volc
      integer ie0,ic1,ic2,ic3,if0,if1,if2,if3,ir1,ir2,ir3
c
c***********************************************************************
c
      volf=dxf(1)*dxf(2)*dxf(3)
      volc=dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=zero
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ic3=ifirstc3,ilastc3
                  if3=ic3*ratio(3)+ir3
                  do ic2=ifirstc2,ilastc2
                     if2=ic2*ratio(2)+ir2
                     do ic1=ifirstc1,ilastc1
                        if1=ic1*ratio(1)+ir1
                        do ie0=ifirstc0,ilastc0+1
                           if0=ie0*ratio(0)
                           arrayc(ie0,ic1,ic2,ic3)=
     &                        arrayc(ie0,ic1,ic2,ic3)
     &                           +arrayf(if0,if1,if2,if3)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=
     &               arrayc(ie0,ic1,ic2,ic3)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot4d1(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      real
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo3:fihi3,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo3:cihi3,
     &          cilo0:cihi0)
      double precision volf,volc
      integer ie1,ic2,ic3,ic0,if1,if2,if3,if0,ir2,ir3,ir0
c
c***********************************************************************
c
      volf=dxf(2)*dxf(3)*dxf(0)
      volc=dxc(2)*dxc(3)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=zero
               enddo
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir3=0,ratio(3)-1
            do ir2=0,ratio(2)-1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ie1=ifirstc1,ilastc1+1
                           if1=ie1*ratio(1)
                           arrayc(ie1,ic2,ic3,ic0)=
     &                        arrayc(ie1,ic2,ic3,ic0)
     &                           +arrayf(if1,if2,if3,if0)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=
     &               arrayc(ie1,ic2,ic3,ic0)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot4d2(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      real
     &  arrayf(filo2:fihi2+1,
     &          filo3:fihi3,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo3:cihi3,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision volf,volc
      integer ie2,ic3,ic0,ic1,if2,if3,if0,if1,ir3,ir0,ir1
c
c***********************************************************************
c
      volf=dxf(3)*dxf(0)*dxf(1)
      volc=dxc(3)*dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=zero
               enddo
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ir3=0,ratio(3)-1
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ic0=ifirstc0,ilastc0
                     if0=ic0*ratio(0)+ir0
                     do ic3=ifirstc3,ilastc3
                        if3=ic3*ratio(3)+ir3
                        do ie2=ifirstc2,ilastc2+1
                           if2=ie2*ratio(2)
                           arrayc(ie2,ic3,ic0,ic1)=
     &                        arrayc(ie2,ic3,ic0,ic1)
     &                           +arrayf(if2,if3,if0,if1)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=
     &               arrayc(ie2,ic3,ic0,ic1)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot4d3(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      real
     &  arrayf(filo3:fihi3+1,
     &          filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo3:cihi3+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision volf,volc
      integer ie3,ic0,ic1,ic2,if3,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      volf=dxf(0)*dxf(1)*dxf(2)
      volc=dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=zero
               enddo
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        do ie3=ifirstc3,ilastc3+1
                           if3=ie3*ratio(3)
                           arrayc(ie3,ic0,ic1,ic2)=
     &                        arrayc(ie3,ic0,ic1,ic2)
     &                           +arrayf(if3,if0,if1,if2)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=
     &               arrayc(ie3,ic0,ic1,ic2)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 4d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx4d0(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2,
     &          filo3:fihi3),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2,
     &          cilo3:cihi3)
      double precision volf,volc
      integer ie0,ic1,ic2,ic3,if0,if1,if2,if3,ir1,ir2,ir3
c
c***********************************************************************
c
      volf=dxf(1)*dxf(2)*dxf(3)
      volc=dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=cmplx(zero,zero)
               enddo
            enddo
         enddo
      enddo

      do ir3=0,ratio(3)-1
         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ic3=ifirstc3,ilastc3
                  if3=ic3*ratio(3)+ir3
                  do ic2=ifirstc2,ilastc2
                     if2=ic2*ratio(2)+ir2
                     do ic1=ifirstc1,ilastc1
                        if1=ic1*ratio(1)+ir1
                        do ie0=ifirstc0,ilastc0+1
                           if0=ie0*ratio(0)
                           arrayc(ie0,ic1,ic2,ic3)=
     &                        arrayc(ie0,ic1,ic2,ic3)
     &                           +arrayf(if0,if1,if2,if3)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ie0=ifirstc0,ilastc0+1
                  arrayc(ie0,ic1,ic2,ic3)=
     &               arrayc(ie0,ic1,ic2,ic3)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx4d1(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double complex
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo3:fihi3,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo3:cihi3,
     &          cilo0:cihi0)
      double precision volf,volc
      integer ie1,ic2,ic3,ic0,if1,if2,if3,if0,ir2,ir3,ir0
c
c***********************************************************************
c
      volf=dxf(2)*dxf(3)*dxf(0)
      volc=dxc(2)*dxc(3)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=cmplx(zero,zero)
               enddo
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir3=0,ratio(3)-1
            do ir2=0,ratio(2)-1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ic3=ifirstc3,ilastc3
                     if3=ic3*ratio(3)+ir3
                     do ic2=ifirstc2,ilastc2
                        if2=ic2*ratio(2)+ir2
                        do ie1=ifirstc1,ilastc1+1
                           if1=ie1*ratio(1)
                           arrayc(ie1,ic2,ic3,ic0)=
     &                        arrayc(ie1,ic2,ic3,ic0)
     &                           +arrayf(if1,if2,if3,if0)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic3=ifirstc3,ilastc3
            do ic2=ifirstc2,ilastc2
               do ie1=ifirstc1,ilastc1+1
                  arrayc(ie1,ic2,ic3,ic0)=
     &               arrayc(ie1,ic2,ic3,ic0)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx4d2(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double complex
     &  arrayf(filo2:fihi2+1,
     &          filo3:fihi3,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo3:cihi3,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision volf,volc
      integer ie2,ic3,ic0,ic1,if2,if3,if0,if1,ir3,ir0,ir1
c
c***********************************************************************
c
      volf=dxf(3)*dxf(0)*dxf(1)
      volc=dxc(3)*dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=cmplx(zero,zero)
               enddo
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ir3=0,ratio(3)-1
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ic0=ifirstc0,ilastc0
                     if0=ic0*ratio(0)+ir0
                     do ic3=ifirstc3,ilastc3
                        if3=ic3*ratio(3)+ir3
                        do ie2=ifirstc2,ilastc2+1
                           if2=ie2*ratio(2)
                           arrayc(ie2,ic3,ic0,ic1)=
     &                        arrayc(ie2,ic3,ic0,ic1)
     &                           +arrayf(if2,if3,if0,if1)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic3=ifirstc3,ilastc3
               do ie2=ifirstc2,ilastc2+1
                  arrayc(ie2,ic3,ic0,ic1)=
     &               arrayc(ie2,ic3,ic0,ic1)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx4d3(
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ifirstc3,
     &  ilastc0,ilastc1,ilastc2,ilastc3,
     &  filo0,filo1,filo2,filo3,fihi0,fihi1,fihi2,fihi3,
     &  cilo0,cilo1,cilo2,cilo3,cihi0,cihi1,cihi2,cihi3
      integer ratio(0:4-1)
      double precision
     &  dxf(0:4-1),
     &  dxc(0:4-1)
      double complex
     &  arrayf(filo3:fihi3+1,
     &          filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo3:cihi3+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision volf,volc
      integer ie3,ic0,ic1,ic2,if3,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      volf=dxf(0)*dxf(1)*dxf(2)
      volc=dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=cmplx(zero,zero)
               enddo
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        do ie3=ifirstc3,ilastc3+1
                           if3=ie3*ratio(3)
                           arrayc(ie3,ic0,ic1,ic2)=
     &                        arrayc(ie3,ic0,ic1,ic2)
     &                           +arrayf(if3,if0,if1,if2)*volf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               do ie3=ifirstc3,ilastc3+1
                  arrayc(ie3,ic0,ic1,ic2)=
     &               arrayc(ie3,ic0,ic1,ic2)/volc
               enddo
            enddo
         enddo
      enddo
c
      return
      end
