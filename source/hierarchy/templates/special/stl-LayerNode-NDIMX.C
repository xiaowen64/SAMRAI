//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/templates/special/stl-LayerNode-NDIMX.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2155 $
// Modified:	$LastChangedDate: 2008-04-28 09:43:00 -0700 (Mon, 28 Apr 2008) $
// Description:	Template instantiation for STL containers of LayerNode.
//

#include <set>

#include "LayerNode.h"
#include "LayerNode.C"

/*
 * This file instantiates any STL classes that may be needed by SAMRAI.
 */

template class std::set<SAMRAI::hier::LayerNode<NDIM> >;
template std::ostream& SAMRAI::hier::operator<< <NDIM>(std::ostream&, SAMRAI::hier::LayerNode<NDIM> const&);

