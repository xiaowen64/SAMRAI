//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/apputils/templates/special/tbox__List-pdat__IndexDataNode-NDIM_appu__CutCell_NDIM_X.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2005/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include <vector>

#include "IndexData.h"
#include "IndexData.C"
#include "CutCell.h"
#include "CutCell.C"
#include "CellGeometry.h"
#include "CellGeometry.C"


template class std::vector< SAMRAI::pdat::IndexDataNode< NDIM,SAMRAI::appu::CutCell<NDIM>,SAMRAI::pdat::CellGeometry<NDIM> > * >;
template class SAMRAI::pdat::IndexDataNode< NDIM,SAMRAI::appu::CutCell<NDIM>,SAMRAI::pdat::CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexIterator< NDIM,SAMRAI::appu::CutCell<NDIM>,SAMRAI::pdat::CellGeometry<NDIM> >;

