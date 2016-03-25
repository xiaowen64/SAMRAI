//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/IntVector.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
