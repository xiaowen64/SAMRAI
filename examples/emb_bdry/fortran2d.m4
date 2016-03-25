define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c***********************************************************************
c  subroutine tagcells
c
c  inputs:
c     ifirst0, ilast0 - patch extents in 0 (X)
c     ifirst1, ilast1 - patch extents in 1 (Y)
c     fgh0,fgh1       - ghosts for flag array
c     tgh0,tgh1       - ghosts for tag array
c     settagval       - value in flag array that sets the tag
c     tag_value       - tag value assigned in tags array
c     flag[nx,ny]     - flag array
c
c  output:
c     tags[nx,ny]     - tag array
c     
c***********************************************************************
      subroutine tagCells(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  fgh0,fgh1,
     &  tgh0,tgh1,
     &  settagval,
     &  tagval,
     &  flag,
     &  tags) 
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1
      integer fgh0,fgh1,tgh0,tgh1
      integer settagval, tagval
c input arrays:
      integer
     &  flag(CELL2dVECG(ifirst,ilast,fgh))

c output arrays:
      integer
     &  tags(CELL2dVECG(ifirst,ilast,tgh))

c
c***********************************************************************
c
      integer i,j,ii,jj

c
c  Loop over elements
c
      do j = ifirst1, ilast1
        do i = ifirst0, ilast0

           if (flag(i,j).eq.settagval) then
              do jj = j-1, j+1
                 do ii = i-1, i+1
                    if ((ii.ge.ifirst0).and.(ii.le.ilast0)) then
                      if ((jj.ge.ifirst1).and.(jj.le.ilast1)) then
                         tags(ii,jj) = tagval
                      endif
                    endif
                 end do
              end do
           endif
        end do
      end do
c
      return
      end
