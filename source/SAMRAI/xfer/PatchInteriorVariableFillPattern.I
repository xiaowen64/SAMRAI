/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class that fills PatchInteriors only
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Return the string name identifier.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
PatchInteriorVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

/*
 *************************************************************************
 *
 * getStencilWidth() throws an error if called.  Only overridding
 * versions of this method in concrete subclasses should be called.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const hier::IntVector&
PatchInteriorVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getZero(d_dim);
}

}
}
