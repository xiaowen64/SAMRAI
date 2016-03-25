//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/toolbox/templates/special/StringSpecial.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	special template file for strings on SGI with CC
//

#include <string>
using namespace std;

#ifdef __KCC_VERSION
#if __KCC_VERSION > 3400   // KCC v3.4
#include <iomanip>
#endif
#endif

#ifdef HAVE_SPECIAL_STRING_OSTREAM_INSTANTIATION
template ostream& std::operator<<(ostream& os, const std::basic_string<char>&);
#endif

/*
 * << operator for complex type requires explicit instantiation with KCC
 * versions 3.4 and higher.  AMW 2/00
 */
#ifdef __KCC_VERSION
#if __KCC_VERSION > 3400   // KCC v3.4
//template std::basic_ostream< char,std::char_traits<char> > &std::operator <<
//     <double, char, std::char_traits<char> >
//     (std::basic_ostream<char,std::char_traits<char> > &,
//      const std::complex<double> &);
//template std::basic_ostream< char,std::char_traits<char> > &std::operator <<
//     <char, std::char_traits<char> >
//     (std::basic_ostream<char,std::char_traits<char> > &, char);
//template std::basic_ostream< char,std::char_traits<char> > &std::operator <<
//     <char, std::char_traits<char> >
//     (std::basic_ostream<char,std::char_traits<char> > &, const char *);
template std::basic_ostream< char,std::char_traits<char> > &std::flush 
     <char, std::char_traits<char> >
     (std::basic_ostream<char,std::char_traits<char> > &);
//template std::basic_string< char,std::char_traits<char>,allocator<char> > 
//     basic_string<char, std::char_traits<char>, std::allocator<char>>
//     (const std::allocator<char> &);
template __kai::omanip_setfill<char> std::setfill<char>(char);
#endif
#endif

#ifdef __DECCXX
#include <sstream>
#include <fstream>
template class std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >;
template class std::basic_string<char, std::char_traits<char>, std::allocator<char> >;
template class std::basic_stringbuf<char, std::char_traits<char> >;
template class std::basic_ostream<char, std::char_traits<char> >;
template class std::basic_ofstream<char, std::char_traits<char> >;
template class std::basic_ios<char, std::char_traits<char> >;
template class std::basic_streambuf<char, std::char_traits<char> >;
template class std::basic_filebuf<char, std::char_traits<char> >;
#endif

#if defined(__GNUC__) && (__GNUC__ == 3 && __GNUC_MINOR__ >= 2)
// Templates for complex functions
template std::basic_string<char, std::char_traits<char>, std::allocator<char> > std::operator+<char, std::char_traits<char>, std::allocator<char> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&);
#endif

