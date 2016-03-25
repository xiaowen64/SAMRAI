//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Array explicit specializations
//

#include "tbox/Array.h"


namespace SAMRAI {
   namespace tbox {

#ifdef HAVE_BOOL
template <> const bool Array<bool>::s_standard_type = true;
#endif

#ifdef HAVE_CHAR
template <> const bool Array<char>::s_standard_type = true;
#endif

#ifdef HAVE_FLOAT
template <> const bool Array<float>::s_standard_type = true;
#endif

template <> const bool Array<double>::s_standard_type = true;
template <> const bool Array<int>::s_standard_type = true;

}
}
