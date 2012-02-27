/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Coarsening algorithm for data transfer between AMR levels
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
void
CoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const boost::shared_ptr<hier::CoarsenOperator>& opcoarsen,
   const boost::shared_ptr<VariableFillPattern>& var_fill_pattern)
{
   registerCoarsen(dst, src, opcoarsen,
      hier::IntVector::getZero(d_dim), var_fill_pattern);
}

/*
 **************************************************************************
 *
 * Reconfigure coarsen schedule to perform operations in this algorithm.
 *
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
CoarsenAlgorithm::checkConsistency(
   const boost::shared_ptr<CoarsenSchedule>& schedule) const
{
   TBOX_ASSERT(schedule);

   return d_coarsen_classes->classesMatch(schedule->getEquivalenceClasses());
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
CoarsenAlgorithm::getDim() const
{
   return d_dim;
}

}
}
