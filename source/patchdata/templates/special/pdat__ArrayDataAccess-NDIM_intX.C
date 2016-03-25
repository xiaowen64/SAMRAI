//
// File:	$URL$
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2275 $
// Modified:	$LastChangedDate: 2008-07-07 13:39:44 -0700 (Mon, 07 Jul 2008) $
// Description:	Special instantiation file for ArrayDataAccess.
//

#include "ArrayDataAccess.h"

namespace SAMRAI {
   namespace pdat {

template MDA_Access<int,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,int>(
   ArrayData<NDIM,int> &array_data,
   int depth );

template MDA_AccessConst<int,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,int>(
   const ArrayData<NDIM,int> &array_data,
   int depth );

}
}

template class MDA_Access<int,NDIM,MDA_OrderColMajor<NDIM> >;
