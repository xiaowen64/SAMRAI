c    
c File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/plotting/appu_cpdat2buf3d.f $
c Package:     SAMRAI application utilities
c Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c Revision:    $LastChangedRevision: 1917 $
c Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c Description: copies data from 3D fortran array to 1D double buffer
c

      subroutine cpfdat2buf3d(
     & gidxlo0, gidxlo1, gidxlo2,
     & bidxlo0, bidxlo1, bidxlo2,
     & bidxhi0, bidxhi1, bidxhi2,
     & gidxhi0, gidxhi1, gidxhi2,
     & farray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1, gidxlo2,
     &        bidxlo0, bidxlo1, bidxlo2,
     &        bidxhi0, bidxhi1, bidxhi2,
     &        gidxhi0, gidxhi1, gidxhi2,
     &        bufsize
      real farray(gidxlo0:gidxhi0,
     &             gidxlo1:gidxhi1,
     &             gidxlo2:gidxhi2)
      double precision buffer(0:bufsize-1)
      integer in0,in1,in2,mark
c     =============================================================

      mark = 0 
      do in2=bidxlo2,bidxhi2
         do in1=bidxlo1,bidxhi1
            do in0=bidxlo0,bidxhi0
            buffer(mark) = dble(farray(in0,in1,in2))
            mark = mark + 1
            enddo
         enddo
      enddo

      return
      end

      subroutine cpddat2buf3d(
     & gidxlo0, gidxlo1, gidxlo2,
     & bidxlo0, bidxlo1, bidxlo2,
     & bidxhi0, bidxhi1, bidxhi2,
     & gidxhi0, gidxhi1, gidxhi2,
     & darray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1, gidxlo2,
     &        bidxlo0, bidxlo1, bidxlo2,
     &        bidxhi0, bidxhi1, bidxhi2,
     &        gidxhi0, gidxhi1, gidxhi2,
     &        bufsize
      double precision darray(gidxlo0:gidxhi0,
     &                        gidxlo1:gidxhi1,
     &                        gidxlo2:gidxhi2)
      double precision buffer(0:bufsize-1)
      integer in0,in1,in2,mark 
c     =============================================================

      mark = 0 
      do in2=bidxlo2,bidxhi2
         do in1=bidxlo1,bidxhi1
            do in0=bidxlo0,bidxhi0
            buffer(mark) = darray(in0,in1,in2)
            mark = mark + 1
            enddo
         enddo
      enddo

      return
      end

      subroutine cpidat2buf3d(
     & gidxlo0, gidxlo1, gidxlo2,
     & bidxlo0, bidxlo1, bidxlo2,
     & bidxhi0, bidxhi1, bidxhi2,
     & gidxhi0, gidxhi1, gidxhi2,
     & iarray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1, gidxlo2,
     &        bidxlo0, bidxlo1, bidxlo2,
     &        bidxhi0, bidxhi1, bidxhi2,
     &        gidxhi0, gidxhi1, gidxhi2,
     &        bufsize,
     &        iarray(gidxlo0:gidxhi0,
     &               gidxlo1:gidxhi1,
     &               gidxlo2:gidxhi2)
      double precision buffer(0:bufsize-1)
      integer in0,in1,in2,mark 
c     =============================================================

      mark = 0 
      do in2=bidxlo2,bidxhi2
         do in1=bidxlo1,bidxhi1
            do in0=bidxlo0,bidxhi0
            buffer(mark) = dble(iarray(in0,in1,in2))
            mark = mark + 1
            enddo
         enddo
      enddo

      return
      end
