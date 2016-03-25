//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/IntVector.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	A n-dimensional integer vector
//

#ifndef included_hier_IntVector_C
#define included_hier_IntVector_C

#include "IntVector.h"

#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "IntVector.I"
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> std::istream& operator >> (std::istream& s, IntVector<DIM>& rhs)
{
   while (s.get() != '(');

   for(int i = 0; i < DIM; i++)
   {
      s >> rhs(i);
      if (i < DIM - 1)
	 while (s.get() != ',');
   }

   while (s.get() != ')');

   return(s); 
}

template<int DIM> std::ostream& operator << (std::ostream& s, 
const IntVector<DIM>& rhs)
{
   s << '(';
   
   for(int i = 0; i < DIM; i++)
   {
      s << rhs(i);
      if (i < DIM - 1)
	 s << ",";
   }
   s << ')';

   return(s);
}

}
}

#endif
