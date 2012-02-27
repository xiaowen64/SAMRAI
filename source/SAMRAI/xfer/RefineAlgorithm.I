/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Refine algorithm for data transfer between AMR levels
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<RefineClasses>&
RefineAlgorithm::getEquivalenceClasses() const
{
   return d_refine_classes;
}

SAMRAI_INLINE_KEYWORD
void
RefineAlgorithm::setEquivalenceClasses(
   const boost::shared_ptr<RefineClasses>& refine_classes)
{
   d_refine_classes = refine_classes;
}

/*
 *************************************************************************
 *
 * Return the dimension
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
RefineAlgorithm::getDim() const
{
   return d_dim;
}

}
}
