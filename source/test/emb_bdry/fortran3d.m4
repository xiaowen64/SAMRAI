define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

c***********************************************************************
c***********************************************************************
c  subroutine tagcells
c
c  inputs:
c     ifirst0, ilast0 - patch extents in 0 (X)
c     ifirst1, ilast1 - patch extents in 1 (Y)
c     ifirst2, ilast2 - patch extents in 2 (Z)
c     fgh0,fgh1,fgh2  - ghosts for flag array
c     tgh0,tgh1,tgh2  - ghosts for tag array
c     settagval       - value in flag array that sets the tag
c     tag_value       - tag value assigned in tags array
c     flag            - flag array
c
c  output:
c     tags            - tag array
c     
c***********************************************************************
      subroutine tagCells(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  fgh0,fgh1,fgh2,
     &  tgh0,tgh1,tgh2,
     &  settagval,
     &  tagval,
     &  flag,
     &  tags) 
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer fgh0,fgh1,fgh2,tgh0,tgh1,tgh2
      integer settagval, tagval
c input arrays:
      integer
     &  flag(CELL3dVECG(ifirst,ilast,fgh))

c output arrays:
      integer
     &  tags(CELL3dVECG(ifirst,ilast,tgh))

c
c***********************************************************************
c
      integer i,j,k,ii,jj,kk

c
c  Loop over elements
c
      do k = ifirst2, ilast2
        do j = ifirst1, ilast1
          do i = ifirst0, ilast0

            if (flag(i,j,k).eq.settagval) then
              do kk = k-1, k+1
                do jj = j-1, j+1
                  do ii = i-1, i+1
                    if ((ii.ge.ifirst0).and.(ii.le.ilast0)) then
                      if ((jj.ge.ifirst1).and.(jj.le.ilast1)) then
                        if ((kk.ge.ifirst2).and.(kk.le.ilast2)) then
                         tags(ii,jj,kk) = tagval
                        endif
                      endif
                    endif
                  end do
                end do
              end do
            endif
          end do
        end do
      end do
c
      return
      end
