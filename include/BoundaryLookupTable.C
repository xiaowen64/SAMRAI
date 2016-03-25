//
// File:        BoundaryLookupTable.C
// Package:     SAMRAI hierarchy 
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:  Lookup table to aid in BoundaryBox construction
//

#ifndef included_hier_BoundaryLookupTable_C
#define included_hier_BoundaryLookupTable_C

#include "BoundaryLookupTable.h"

#include "tbox/ShutdownRegistry.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#include "BoundaryLookupTable.I"
#endif


namespace SAMRAI {
   namespace hier {


template<int DIM> BoundaryLookupTable<DIM>*
BoundaryLookupTable<DIM>::s_lookup_table_instance =
   (BoundaryLookupTable<DIM>*) NULL;
template<int DIM> bool BoundaryLookupTable<DIM>::s_registered_callback = false;
template<int DIM> bool BoundaryLookupTable<DIM>::s_using_original_locations = false;


/*
*************************************************************************
*                                                                       *
* Lookup table constructor and destructor.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>  BoundaryLookupTable<DIM>::BoundaryLookupTable()
{

   /*
    * These "fake" uses of statics avoid the Intel compiler
    * optimizing them away.  They are only used in 
    * inline methods.
    */
   if(s_lookup_table_instance == NULL) {
      s_lookup_table_instance = NULL;
   }
   if(s_registered_callback == false) {
      s_lookup_table_instance = false;
   }

   if (d_table[0].isNull()) {
      int factrl[DIM+1] = {1};
      for (int i = 1; i <= DIM; i++) factrl[i] = i * factrl[i-1];
      d_ncomb.resizeArray(DIM);
      d_max_li.resizeArray(DIM);
      for (int codim = 1; codim <= DIM; codim++) {
	 int cdm1 = codim-1;
         d_ncomb[cdm1] = factrl[DIM]/(factrl[codim]*factrl[DIM-codim]);

	 tbox::Array<int> work;
	 work.resizeArray(codim*d_ncomb[cdm1]);
	 buildTable(work.getPointer(), codim, 1);

	 d_table[cdm1].resizeArray(d_ncomb[cdm1]);
	 for (int j=0; j<d_ncomb[cdm1]; j++) {
	    d_table[cdm1][j].resizeArray(codim);
	    for (int k=0; k<codim; k++) {
	       d_table[cdm1][j][k] = work[j*codim+k]-1;
	    }
	 }

         d_max_li[cdm1] = d_ncomb[cdm1]*(1<<codim);
      }
   }
}

template<int DIM>  BoundaryLookupTable<DIM>::~BoundaryLookupTable()
{
}

/*
*************************************************************************
*                                                                       *
* Recursive function that computes the combinations in the lookup       *
* table.                                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoundaryLookupTable<DIM>::buildTable(int *table, int codim, int ibeg)
{
   static int work[DIM];
   static int lvl = 0;
   static int *ptr;
   lvl++;                                  
   if (lvl == 1) ptr = table;              
   int iend = DIM - codim + lvl;          
   for (int i = ibeg; i <= iend; i++) {    
      work[lvl-1] = i;
      if (lvl != codim) {
	 buildTable(ptr, codim, i+1);
      } else {
         for (int j = 0; j < codim; j++) {
	    *(ptr+j) = work[j];
	 }
	 ptr += codim;
      }
   }
   lvl--;
}

template<int DIM> void BoundaryLookupTable<DIM>::setUsingOriginalLocations(
   const bool use_original)
{
   s_using_original_locations = use_original;
}
 
}
}

#endif
