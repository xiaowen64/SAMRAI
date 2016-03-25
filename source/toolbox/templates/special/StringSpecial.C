//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/templates/special/StringSpecial.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1959 $
// Modified:	$LastChangedDate: 2008-02-06 13:52:33 -0800 (Wed, 06 Feb 2008) $
// Description:	special template file for strings on SGI with CC
//

#include <string>
using namespace std;

#ifdef HAVE_SPECIAL_STRING_OSTREAM_INSTANTIATION
template ostream& std::operator<<(ostream& os, const std::basic_string<char>&);
#endif

