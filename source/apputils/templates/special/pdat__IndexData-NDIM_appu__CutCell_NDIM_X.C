//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/templates/special/pdat__IndexData-NDIM_appu__CutCell_NDIM_X.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Special class for index data
//

#include "IndexData.h"
#include "IndexData.C"
#include "CutCell.h"

namespace SAMRAI {
   namespace pdat {
template class pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> >;
template class pdat::IndexIterator< NDIM,appu::CutCell<NDIM> >;
}
}
