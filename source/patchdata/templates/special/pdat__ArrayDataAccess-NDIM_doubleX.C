//
// File:	$URL$
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2244 $
// Modified:	$LastChangedDate: 2008-07-02 11:10:50 -0700 (Wed, 02 Jul 2008) $
// Description:	Special instantiation file for ArrayDataAccess.
//

#include "ArrayDataAccess.h"

namespace SAMRAI {
   namespace pdat {

template MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,double>(
   ArrayData<NDIM,double> &array_data,
   int depth );

template MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> >
ArrayDataAccess::access<NDIM,double>(
   const ArrayData<NDIM,double> &array_data,
   int depth );

}
}

template class MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >;
