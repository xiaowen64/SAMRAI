//
// File:	$Id: stl-LayerNode-NDIMX.C 346 2005-05-09 19:43:12Z gunney $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 346 $
// Modified:	$Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
// Description:	Template instantiation for STL containers of LayerNode.
//

#include <set>
#include <vector>
#include "LayerNode.h"
#include "LayerNode.C"

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

template class std::set<SAMRAI::hier::LayerNode<NDIM> >;

#if defined(__xlC__)

template std::ostream& SAMRAI::hier::operator<< <NDIM>(std::ostream&, SAMRAI::hier::LayerNode<NDIM> const&);


#elif defined(__INTEL_COMPILER)

template std::ostream& SAMRAI::hier::operator<< <NDIM>(std::ostream&, SAMRAI::hier::LayerNode<NDIM> const&);


#elif defined(__GNUC__)


template class std::_Rb_tree<int, std::pair<int const, std::set<SAMRAI::hier::LayerNode<NDIM>, std::less<SAMRAI::hier::LayerNode<NDIM> >, std::allocator<SAMRAI::hier::LayerNode<NDIM> > > >, std::_Select1st<std::pair<int const, std::set<SAMRAI::hier::LayerNode<NDIM>, std::less<SAMRAI::hier::LayerNode<NDIM> >, std::allocator<SAMRAI::hier::LayerNode<NDIM> > > > >, std::less<int>, std::allocator<std::pair<int const, std::set<SAMRAI::hier::LayerNode<NDIM>, std::less<SAMRAI::hier::LayerNode<NDIM> >, std::allocator<SAMRAI::hier::LayerNode<NDIM> > > > > >;
template std::basic_ostream<char, std::char_traits<char> >& SAMRAI::hier::operator<< <NDIM>(std::basic_ostream<char, std::char_traits<char> >&, SAMRAI::hier::LayerNode<NDIM> const&);
template class std::_Rb_tree<SAMRAI::hier::LayerNode<NDIM>, SAMRAI::hier::LayerNode<NDIM>, std::_Identity<SAMRAI::hier::LayerNode<NDIM> >, std::less<SAMRAI::hier::LayerNode<NDIM> >, std::allocator<SAMRAI::hier::LayerNode<NDIM> > >;

#endif
