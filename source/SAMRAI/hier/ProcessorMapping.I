/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   tbox 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
ProcessorMapping::~ProcessorMapping()
{
}

SAMRAI_INLINE_KEYWORD
void ProcessorMapping::setNumberNodes(
   const int n)
{
   d_nodes = n;
}

SAMRAI_INLINE_KEYWORD
int ProcessorMapping::getProcessorAssignment(
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));

   return d_mapping[i];
}

SAMRAI_INLINE_KEYWORD
void ProcessorMapping::setProcessorAssignment(
   const int i,
   const int p)
{
   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
   TBOX_ASSERT((p >= 0) && (p < d_nodes));

   d_mapping[i] = p % d_nodes;
}

SAMRAI_INLINE_KEYWORD
tbox::Array<int> ProcessorMapping::getProcessorMapping() const
{
   return d_mapping;
}

SAMRAI_INLINE_KEYWORD
int ProcessorMapping::getSizeOfMappingArray() const
{
   return d_mapping.getSize();
}

SAMRAI_INLINE_KEYWORD
bool ProcessorMapping::isMappingLocal(
   const int i) const
{

   TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));

   return d_mapping[i] == d_my_rank;
}

}
}
