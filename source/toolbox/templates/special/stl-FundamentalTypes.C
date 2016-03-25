//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/templates/special/stl-FundamentalTypes.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2150 $
// Modified:	$LastChangedDate: 2008-04-24 13:44:21 -0700 (Thu, 24 Apr 2008) $
// Description:	Template instantiation for STL containers of int.
//

#include <vector>
#include <set>
#include <map>

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


template class std::set<int>;


template class std::map<int,int>;


template class std::map<int,std::vector<int> >;

/*
 * IBM XLC fails on this.
 */
#if !defined(__xlC__)
template void std::set<int>::insert<int*>(int *, int *);
template void std::set<int>::insert<const int *>(const int *, const int *);
#endif


template class std::vector<char>;
template void std::vector<char>::insert<char*>(std::vector<char>::iterator, char *, char *);
template void std::vector<char>::insert<const char *>(std::vector<char>::iterator, const char *, const char *);

template class std::vector<int>;
template void std::vector<int>::insert<int*>(std::vector<int>::iterator, int *, int *);
template void std::vector<int>::insert<const int *>(std::vector<int>::iterator, const int *, const int *);

template class std::vector<float>;
template void std::vector<float>::insert<float*>(std::vector<float>::iterator, float *, float *);
template void std::vector<float>::insert<const float *>(std::vector<float>::iterator, const float *, const float *);

template class std::vector<double>;
template void std::vector<double>::insert<double*>(std::vector<double>::iterator, double *, double *);
template void std::vector<double>::insert<const double *>(std::vector<double>::iterator, const double *, const double *);




