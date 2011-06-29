/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:    
 *
 ************************************************************************/

#include "SAMRAI/pdat/ArrayDataAccess.h"

template<int DIM>
MDA_Access<double, DIM, MDA_OrderColMajor<DIM> > access(
   SAMRAI::pdat::ArrayData<double>& array_data,
   int depth)
{
   MDA_Access<double, DIM, MDA_OrderColMajor<DIM> >
   r(array_data.getPointer(depth),
     (const int *)array_data.getBox().lower(),
     (const int *)array_data.getBox().upper());
   return r;
}

template<int DIM>
const MDA_AccessConst<double, DIM, MDA_OrderColMajor<DIM> > access(
   const SAMRAI::pdat::ArrayData<double>& array_data,
   int depth)
{
   MDA_AccessConst<double, DIM, MDA_OrderColMajor<DIM> >
   r(array_data.getPointer(depth),
     (const int *)array_data.getBox().lower(),
     (const int *)array_data.getBox().upper());
   return r;
}
