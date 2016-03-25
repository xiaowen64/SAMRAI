//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/templates/special/ComplexSpecial.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2004 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1957 $
// Modified:	$LastChangedDate: 2008-02-05 16:31:08 -0800 (Tue, 05 Feb 2008) $
// Description:	special template file for complex numbers on SGI with CC
//

#include "tbox/Complex.h"

#ifdef HAVE_SPECIAL_COMPLEX_OSTREAM_INSTANTIATION
template ostream& std::operator<<(ostream&,const std::complex<double>&);
#endif

// Templates for complex functions, these are not always required but 
// should always be valid.
template double std::abs<double>(std::complex<double> const&);
template std::complex<double> std::polar<double>(double const&, double const&);
template std::complex<double> std::conj<double>(std::complex<double> const&);
template std::complex<double> std::cos<double>(std::complex<double> const&);
template std::complex<double> std::cosh<double>(std::complex<double> const&);
template std::complex<double> std::exp<double>(std::complex<double> const&);
template std::complex<double> std::log<double>(std::complex<double> const&);
template std::complex<double> std::sin<double>(std::complex<double> const&);
template std::complex<double> std::sinh<double>(std::complex<double> const&);
template std::complex<double> std::sqrt<double>(std::complex<double> const&);
template double std::norm<double>(std::complex<double> const&);

