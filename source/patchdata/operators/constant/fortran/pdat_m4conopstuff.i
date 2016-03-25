c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/operators/constant/fortran/pdat_m4conopstuff.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
