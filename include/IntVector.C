//
// File:	IntVector.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

template<int DIM> istream& operator >> (istream& s, IntVector<DIM>& rhs)
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

template<int DIM> ostream& operator << (ostream& s, 
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
