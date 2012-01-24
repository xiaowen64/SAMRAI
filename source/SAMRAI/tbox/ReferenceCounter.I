/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Reference counting class for Array
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
ReferenceCounter::ReferenceCounter()
{
   d_references = 1;
   d_next = (ReferenceCounter *)NULL;
}

SAMRAI_INLINE_KEYWORD
ReferenceCounter::~ReferenceCounter()
{
   if ((d_next) && (--d_next->d_references == 0)) {
      delete d_next;
   }
}

SAMRAI_INLINE_KEYWORD
bool ReferenceCounter::deleteReference()
{
   return --d_references == 0;
}

SAMRAI_INLINE_KEYWORD
void ReferenceCounter::addReference()
{
   d_references++;
}

}
}
