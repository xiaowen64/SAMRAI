/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Auto-release handle to prevent using invalid Connector data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
BoxLevelHandle::isAttached() const
{
   return d_box_level != NULL;
}

/*
 ************************************************************************
 * To be called by the object's BoxLevel when the
 * BoxLevel changes in a way that can invalidate the
 * Connector data.
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelHandle::detachMyBoxLevel()
{
   d_box_level = NULL;
}

}
}
