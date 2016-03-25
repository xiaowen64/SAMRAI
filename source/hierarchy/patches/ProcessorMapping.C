//
// File:   ProcessorMapping.C
// Package:   SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:   $Revision: 173 $
// Modified:   $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:   tbox::Array of processor mappings of patches to processors
//

#include "ProcessorMapping.h"
#include "tbox/MPI.h"

#ifdef DEBUG_NO_INLINE
#include "ProcessorMapping.I"
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace hier {

ProcessorMapping::ProcessorMapping()
:  d_my_rank(tbox::MPI::getRank()),
   d_nodes(tbox::MPI::getNodes()),
   d_mapping(0),
   d_local_index_count(-1)
{
}

ProcessorMapping::ProcessorMapping(const int n)
:  d_my_rank(tbox::MPI::getRank()),
   d_nodes(tbox::MPI::getNodes()),
   d_mapping(n),
   d_local_index_count(-1)
{
   for (int i = 0; i < n; i++) {
      d_mapping[i] = 0;
   }
}

ProcessorMapping::ProcessorMapping(
   const ProcessorMapping& mapping)
:  d_my_rank(tbox::MPI::getRank()),
   d_nodes(tbox::MPI::getNodes()),
   d_mapping(mapping.d_mapping.getSize()),
   d_local_index_count(-1)
{
   const int n = d_mapping.getSize();
   for (int i = 0; i < n; i++) {
      d_mapping[i] = mapping.d_mapping[i];
   }
}

ProcessorMapping::ProcessorMapping(
   const tbox::Array<int>& mapping)
:  d_my_rank(tbox::MPI::getRank()),
   d_nodes(tbox::MPI::getNodes()),
   d_local_index_count(-1)
{
   setProcessorMapping(mapping);
}

void ProcessorMapping::setMappingSize(const int n) 
{
   d_mapping.resizeArray(n);

   for (int i = 0; i < n; i++) {
      d_mapping[i] = 0;
   }
   d_local_index_count = -1;
}

void ProcessorMapping::setProcessorMapping(const tbox::Array<int>& mapping)
{
   d_mapping.resizeArray(mapping.getSize());

   for (int i = 0; i < d_mapping.getSize(); i++) {
      //  (mapping[i] % d_nodes) keeps patches from being assigned
      //  non-existent processors.
      setProcessorAssignment(i,mapping[i] % d_nodes); 
   }
   d_local_index_count = -1;
}

int ProcessorMapping::getNumberOfLocalIndices() const
{
   computeLocalIndices();
   return( d_local_index_count );
}

const tbox::Array<int>& ProcessorMapping::getLocalIndices() const
{
   computeLocalIndices();
   return(d_local_indices);
}

void ProcessorMapping::computeLocalIndices() const
{
   if (d_local_index_count != -1) {
      return;
   }

   /*
    * first, count the number of local indices,
    * so we can set the array size.
    */
   const int n = d_mapping.getSize();
   d_local_index_count = 0;

   for (int i = 0; i < n; i++) {
      if (d_mapping[i] == d_my_rank) {
         ++d_local_index_count;
      }
   }

   /*
    * second, resize the array and fill in the data 
    */
   d_local_indices.resizeArray(d_local_index_count);
   int idx = 0;
   for (int i = 0; i < n; i++) {
      if (d_mapping[i] == d_my_rank) {
         d_local_indices[idx++] = i;
      }
   }
}

}
}
