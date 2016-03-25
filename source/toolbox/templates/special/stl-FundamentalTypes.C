//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/templates/special/stl-FundamentalTypes.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Template instantiation for STL containers of int.
//

#include <vector>
#include <set>

/*
 * This file instantiates STL classes and functions
 * for template arguments that are just fundamental
 * types (int, double, etc.).
 *
 * This comment on adding explicit instantiations
 * also applies for instantiating STL for template
 * arguments that are classes and structs.
 *
 * Motivation:
 *
 * To avoid long link times, the SAMRAI configuration
 * disables implicit template instantation.  Instantiations
 * are done explicitly so that each unique one is done
 * just once.
 *
 * A potential difficulty with explicitly instantiating
 * STL is that *all* templates used must be instantiated
 * even if they are used internally (i.e., not in the
 * public interface of the STL).  It is impossible to
 * know a priori the internally used templates, as this
 * is implementationally dependent.
 *
 * We determine what must be explicitly instantiated
 * by seeing what is unresolved at link time.
 * Unfortunately, we know of no other method for doing this!
 *
 * How to add STL template instantiations:
 *
 * When an STL template or a template used internally by STL
 * comes up missing at link time, it should be added to one
 * of the source/<package>/templates/special/stl-*.C files
 * in SAMRAI.  Examples:
 * - Symbols resulting from instantiating STL with fundamental
 *   types (int, double, etc) are explicitly instantiated in
 *   "source/toolbox/templates/special/stl-FundamentalTypes.C"
 * - Symbols resulting from instantiating STL with LayerNode<NDIM>
 *   are explicitly instantiated in
 *   "source/mesh/templates/special/stl-LayerNode-NDIMX.C"
 *
 * Because the specific missing symbols may be implementation-
 * dependent, these files have #if directives that
 * activate just the specific code required by specific systems.
 * Important note: The Intel compiler defines __GNUC__ by default.
 * Therefore it is safest to set up the #if blocks to handle each
 * system exclusively and check the Intel compiler case before
 * checking the GNU case.
 *
 * To add explicit instantiation code, it helps if you understand
 * the syntax required.  Instantiating a class seems to instantiate
 * all non-templated member functions, so non-templated member
 * functions need not be individually instantiated.  Templated
 * functions and member functions must be instantiated individually.
 */

template class std::vector<int>;

template class std::set<int>;

#if defined(__INTEL_COMPILER)

#elif defined(__GNUC__)

template void std::fill<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int>(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int const&);

template int* std::fill_n<int*, unsigned, int>(int*, unsigned, int const&);

template class __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > > std::fill_n<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, unsigned, int>(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, unsigned, int const&);

template class std::_Rb_tree<int, int, std::_Identity<int>, std::less<int>, std::allocator<int> >;

template class std::_Rb_tree<int, std::pair<int const, std::vector<int, std::allocator<int> > >, std::_Select1st<std::pair<int const, std::vector<int, std::allocator<int> > > >, std::less<int>, std::allocator<std::pair<int const, std::vector<int, std::allocator<int> > > > >;

// Symbol depends on version of GNU compiler that is being used
#if __GNUC__ > 3 
   template void std::_Rb_tree<int, int, std::_Identity<int>, std::less<int>, std::allocator<int> >::_M_insert_unique<int*>(int*, int*);
   template void std::fill<int*, int>(int*, int*, int const&);
#else
   template void std::_Rb_tree<int, int, std::_Identity<int>, std::less<int>, std::allocator<int> >::insert_unique<int*>(int*, int*);
#endif


// The following are required by the gps-gcc test case.

template
__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > > std::fill_n<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, unsigned long, int>(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, unsigned long, int const&);

template
int* std::fill_n<int*, unsigned long, int>(int*, unsigned long, int const&);

// The following is required when configuring for shared SAMRAI libraries.

template void std::vector<int, std::allocator<int> >::_M_range_insert<int*>(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int*, int*, std::forward_iterator_tag);

#endif
