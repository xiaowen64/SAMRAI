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

template MDA_Access<char,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,char>(
   ArrayData<NDIM,char> &array_data,
   int depth );

template MDA_AccessConst<char,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,char>(
   const ArrayData<NDIM,char> &array_data,
   int depth );

}
}

template class MDA_Access<char,NDIM,MDA_OrderColMajor<NDIM> >;
