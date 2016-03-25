//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/toolbox/templates/special/StringInst.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	special template file for strings on SGI with CC
//

#include <string>
using namespace std;

#ifdef HAVE_SPECIAL_STRING_OSTREAM_INSTANTIATION
template ostream& std::operator<<(ostream& os, const std::basic_string<char>&);
#endif

