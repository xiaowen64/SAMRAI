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
bool ReferenceCounter::deleteReference()
{
   return --d_references == 0;
}

SAMRAI_INLINE_KEYWORD
void ReferenceCounter::addReference()
{
   d_references++;
}

SAMRAI_INLINE_KEYWORD
bool ReferenceCounter::isFinalized() {
   return s_is_finalized;
}

}
}
