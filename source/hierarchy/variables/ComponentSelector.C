//
// File:	ComponentSelector.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 407 $
// Modified:	$Date: 2005-06-01 09:49:54 -0700 (Wed, 01 Jun 2005) $
// Description:	Simple bit vector of a fixed length (128 bits)
//

#ifndef included_hier_ComponentSelector_C
#define included_hier_ComponentSelector_C

#include "ComponentSelector.h"

#ifdef DEBUG_NO_INLINE
#include "ComponentSelector.I"
#endif

namespace SAMRAI {
   namespace hier {

int ComponentSelector::s_bits_per_long = 8*sizeof(unsigned long);

void ComponentSelector::printClassData( ostream &os ) const
{
   int i;
   const int number_of_bits = getSize();
   for ( i=0; i<number_of_bits; ++i ) {
      os << " | Bit " << i << " = " << isSet(i);
   }
   os << "|\n";
   return;
}

}
}

#endif
