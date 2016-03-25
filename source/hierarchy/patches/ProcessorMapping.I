//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/ProcessorMapping.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	tbox::Array of processor mappings of patches to processors
//


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
   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
#endif
   return(d_mapping[i]);
}

inline
void ProcessorMapping::setProcessorAssignment(const int i, const int p)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
   TBOX_ASSERT((p >= 0) && (p < d_nodes));
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
   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
#endif
   return(d_mapping[i] == d_my_rank);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

