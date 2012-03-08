/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface to application-specific patch functions in support
 *                Method of Lines integration algorithm
 *
 ************************************************************************/

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
MethodOfLinesPatchStrategy::getInteriorWithGhostsContext() const
{
   return d_interior_with_ghosts;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
MethodOfLinesPatchStrategy::getInteriorContext() const
{
   return d_interior;
}

SAMRAI_INLINE_KEYWORD
void
MethodOfLinesPatchStrategy::setInteriorWithGhostsContext(
   const boost::shared_ptr<hier::VariableContext>& context)
{
   d_interior_with_ghosts = context;
}

SAMRAI_INLINE_KEYWORD
void
MethodOfLinesPatchStrategy::setInteriorContext(
   const boost::shared_ptr<hier::VariableContext>& context)
{
   d_interior = context;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
MethodOfLinesPatchStrategy::getDim() const
{
   return d_dim;
}

}
}
