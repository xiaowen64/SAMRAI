//
// File:	ComplexSpecial.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2004 The Regents of the University of California
// Revision:	$Revision: 488 $
// Modified:	$Date: 2005-07-25 18:15:19 -0700 (Mon, 25 Jul 2005) $
// Description:	special template file for complex numbers on SGI with CC
//

#include "tbox/Complex.h"

#ifdef HAVE_SPECIAL_COMPLEX_OSTREAM_INSTANTIATION
template ostream& std::operator<<(ostream&,const std::complex<double>&);
#endif

/*
 * << operator for complex type requires explicit instantiation with KCC
 * versions 3.4 and higher. 
 */
#if __KCC_VERSION > 3400   
template std::basic_ostream< char,std::char_traits<char> > &std::operator <<
     <double, char, std::char_traits<char> >
     (std::basic_ostream<char,std::char_traits<char> > &,
      const std::complex<double> &);
#endif

#if defined(__SUNPRO_CC) 
// On Sun we need to explicitaly instantiate this template function
template complex<double> conj(const complex<double>&);
#endif

#if defined(__DECCXX)
template complex<double> std::conj(const complex<double>&);
#endif

#if defined(__GNUC__)  && (__GNUC__ == 3 && __GNUC_MINOR__ >= 2)
// Templates for complex functions
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
#if !defined(__INTEL_COMPILER)
template double std::__cmath_power<double>(double, unsigned);
//template std::complex<double> std::pow<double>(std::complex<double> const&, std::complex<double> const&);
template double std::arg<double>(std::complex<double> const&);
#ifdef HAVE_ISNAN_TEMPLATE
template int __gnu_cxx::isnan<float>(float);
template int __gnu_cxx::__capture_isnan<float>(float);
#endif
#endif
#endif
