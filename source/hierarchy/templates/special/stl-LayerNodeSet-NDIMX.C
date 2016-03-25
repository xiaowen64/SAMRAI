//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/templates/special/stl-LayerNodeSet-NDIMX.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Template instantiation for STL containers of LayerNodeSet.
//

#include <set>
#include <vector>
#include "LayerEdgeSet.h"
#include "LayerEdgeSet.C"

/*
 * This file instantiates any STL classes that may be needed by SAMRAI.
 * We determine what must be explicitly instantiated by seeing what is
 * unresolved at link time.  Unfortunately, we know of no other method
 * for doing this!
 *
 * The #if statements activate code required by specific systems.
 *
 * Important note: The Intel compiler defines __GNUC__ by default.
 * Therefore it is safest to set up the #if blocks to handle each
 * system exclusively and check the Intel compiler case before
 * checking the GNU case.
 */

template class std::map<SAMRAI::hier::LayerEdgeSet<NDIM>::LocalIndex,
                        SAMRAI::hier::LayerEdgeSet<NDIM>::NabrContainer >;

#if defined(__xlC__)


#elif defined(__INTEL_COMPILER)


#elif defined(__GNUC__)


template class std::_Rb_tree<int, std::pair<int const, SAMRAI::hier::LayerEdgeSet<NDIM>::CommunicationStruct>, std::_Select1st<std::pair<int const, SAMRAI::hier::LayerEdgeSet<NDIM>::CommunicationStruct> >, std::less<int>, std::allocator<std::pair<int const, SAMRAI::hier::LayerEdgeSet<NDIM>::CommunicationStruct> > >;


#endif
