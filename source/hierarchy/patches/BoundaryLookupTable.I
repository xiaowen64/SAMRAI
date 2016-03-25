//
// File:	BoundaryLookupTable.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 699 $
// Modified:	$Date: 2005-11-03 12:30:16 -0800 (Thu, 03 Nov 2005) $
// Description:	
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "tbox/ShutdownRegistry.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


/*
*************************************************************************
*                                                                       *
* Static lookup table member functions.                                 *
*                                                                       *
*************************************************************************
*/
                                                                                
template<int DIM> inline
BoundaryLookupTable<DIM>* BoundaryLookupTable<DIM>::getLookupTable()
{
   if (!s_lookup_table_instance) {
      s_lookup_table_instance = new BoundaryLookupTable<DIM>();
   }
   if (!s_registered_callback) {
      tbox::ShutdownRegistry::registerShutdownRoutine(freeLookupTable,
	tbox::ShutdownRegistry::priorityBoundaryLookupTable);
      s_registered_callback = true;
   }
   return(s_lookup_table_instance);
}
                                                                                
template<int DIM> inline
void BoundaryLookupTable<DIM>::freeLookupTable()
{
   if (s_lookup_table_instance) delete s_lookup_table_instance;
   s_lookup_table_instance = ((BoundaryLookupTable<DIM>*) NULL);
}


/*
*************************************************************************
*                                                                       *
* Inline accessory functions to retrieve data members.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
const tbox::Array<int>& BoundaryLookupTable<DIM>::getMaxLocationIndices() const
{
   return(d_max_li);
}

template<int DIM> inline
const tbox::Array<int>& BoundaryLookupTable<DIM>::getDirections(
   const int loc,
   const int codim) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((codim > 0) && (codim <= DIM));
   assert((loc >= 0) && (loc < d_max_li[codim-1]));
#endif

   int location_index = loc;
   if (codim == 2) {
      location_index = mapLocationIndex(loc);
   }

   int iloc = location_index/(1<<codim);
   return(d_table[codim-1][iloc]);
}

template<int DIM> inline
bool BoundaryLookupTable<DIM>::isLower(const int loc,
                                        const int codim,
                                        const int index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((codim > 0) && (codim <= DIM));
   assert((loc >= 0) && (loc < d_max_li[codim-1]));
   assert((index >= 0) && (index < codim));
#endif
   return(!isUpper(loc, codim, index));
}

template<int DIM> inline
bool BoundaryLookupTable<DIM>::isUpper(const int loc,
                                        const int codim,
                                        const int index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((codim > 0) && (codim <= DIM));
   assert((loc >= 0) && (loc < d_max_li[codim-1]));
   assert((index >= 0) && (index < codim));
#endif

   int location_index = loc;
   if (codim == 2) {
      location_index = mapLocationIndex(loc);
   }

   return((location_index%(1<<codim)) & (1<<(index)));
}

/*
*************************************************************************
*                                                                       *
* For backward compatibility, map location index from previous system   *
* to a value that is understood by the lookup table.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
int BoundaryLookupTable<DIM>::mapLocationIndex(const int loc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(loc < d_max_li[1]);
   assert(loc >= 0);
#endif

   int retval = loc;
   if (s_using_original_locations && (DIM == 3)) {
                                                                                
      if (loc < 4) {
         retval = loc + 8;
      } else if (loc >= 8) {
         retval = loc - 8;
      } else if ((loc == 4) || (loc == 7)) {
         retval = loc;
      } else {
         retval = (loc == 5) ? 6 : 5;
      }
   }                                                                              
   return(retval);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
