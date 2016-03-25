//
// File:	$Id: stl-LayerNodeSet-NDIMX.C 346 2005-05-09 19:43:12Z gunney $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 346 $
// Modified:	$Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
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
