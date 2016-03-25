c
c  File:        pdat_m4conopstuff.i
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2005 The Regents of the University of California
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: m4 include file for constant patchdata transfer routines.
c
define(coarsen_index,`dnl
         if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl
define(coarsen_face_index,`dnl
         it=2*$1+$3
         if (it.le.0) then
            $2=it/(2*$3)-1
         else
            $2=(it-1)/(2*$3)
         endif
')dnl
