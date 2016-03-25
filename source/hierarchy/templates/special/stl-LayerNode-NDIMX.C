//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/templates/special/stl-LayerNode-NDIMX.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
