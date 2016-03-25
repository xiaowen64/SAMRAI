//
// File:	ProcessorMapping.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	tbox::Array of processor mappings of patches to processors
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

inline
ProcessorMapping::~ProcessorMapping()
{
}

inline
void ProcessorMapping::setNumberNodes(const int n)
{
   d_nodes = n;
}

inline
int ProcessorMapping::getProcessorAssignment(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_mapping.getSize()));
#endif
   return(d_mapping[i]);
}

inline
void ProcessorMapping::setProcessorAssignment(const int i, const int p)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_mapping.getSize()));
   assert((p >= 0) && (p < d_nodes));
#endif
   d_mapping[i] = p % d_nodes;
}

inline
tbox::Array<int> ProcessorMapping::getProcessorMapping() const
{
   return d_mapping;
}

inline
int ProcessorMapping::getSizeOfMappingArray() const
{
   return(d_mapping.getSize());
}

inline
bool ProcessorMapping::isMappingLocal(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_mapping.getSize()));
#endif
   return(d_mapping[i] == d_my_rank);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

