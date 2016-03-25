//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/ComponentSelector.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

void ComponentSelector::printClassData( std::ostream &os ) const
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
