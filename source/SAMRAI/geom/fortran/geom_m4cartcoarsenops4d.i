c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 4d Cartesian coarsen operators
c
define(NDIM,4)dnl
define(SAMRAI_FORTDIR,../../pdat/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim4d.i)dnl
define(cart_coarsen_op_subroutine_head_4d,`dnl
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
      integer ratio(0:NDIM-1)
      double precision
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1)
')dnl
c
define(cart_wgtavg_cell_body_4d,`dnl
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)*dxf(3)
      dVc = dxc(0)*dxc(1)*dxc(2)*dxc(3)

      do ic3=ifirstc3,ilastc3
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
ifelse($1,`double complex',`dnl
                  arrayc(ic0,ic1,ic2,ic3)=cmplx(zero,zero)
',`dnl
                  arrayc(ic0,ic1,ic2,ic3)=zero
')dnl
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
')dnl
c
define(cart_wgtavg_op_cell_4d,`dnl
cart_coarsen_op_subroutine_head_4d()dnl
      $1
     &  arrayf(CELL4d(filo,fihi,0)),
     &  arrayc(CELL4d(cilo,cihi,0))
      double precision dVf,dVc
      integer ic0,ic1,ic2,ic3,if0,if1,if2,if3,ir0,ir1,ir2,ir3
cart_wgtavg_cell_body_4d($1)dnl
')dnl
