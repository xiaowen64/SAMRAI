//
// File:	$URL: svn+ssh://tux262.llnl.gov/usr/casc/samrai/repository/SAMRAI/trunk/source/hierarchy/templates/special/stl-LayerNodeSet-NDIMX.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2155 $
// Modified:	$LastChangedDate: 2008-04-28 09:43:00 -0700 (Mon, 28 Apr 2008) $
// Description:	Template instantiation for STL containers of LayerNodeSet.
//

#include "LayerEdgeSet.h"
#include <map>

/*
 * This file instantiates any STL classes that may be needed by SAMRAI.
 */

template class std::map<SAMRAI::hier::LayerEdgeSet<NDIM>::LocalIndex,
                        SAMRAI::hier::LayerEdgeSet<NDIM>::NabrContainer >;

template class std::map< int, SAMRAI::hier::LayerEdgeSet<NDIM>::CommunicationStruct >;
