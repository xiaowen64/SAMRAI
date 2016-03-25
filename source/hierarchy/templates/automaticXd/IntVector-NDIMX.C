//
// File:	IntVector-NDIMX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "IntVector.h"
#include "IntVector.C"

namespace SAMRAI {
   namespace hier {
template class IntVector< NDIM >;
template istream& operator >> (istream& s, IntVector<NDIM>& box);
template ostream& operator << (ostream& s, const IntVector<NDIM>& box);
}
}
