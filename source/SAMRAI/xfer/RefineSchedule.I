/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Refine schedule for data transfer between AMR levels
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 **************************************************************************
 *
 * Return const pointer to equivalence classes used in schedule.
 *
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<RefineClasses>&
RefineSchedule::getEquivalenceClasses() const
{
   return d_refine_classes;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector&
RefineSchedule::getBoundaryFillGhostWidth() const
{
   return d_boundary_fill_ghost_width;
}

/*
 **************************************************************************
 *
 * Fill the component selector argument with the components needed to
 * allocate source data.
 *
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineSchedule::initializeSourceVector(
   hier::ComponentSelector& allocate_vector) const
{
   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      allocate_vector.setFlag(d_refine_items[iri]->d_src);
   }
}

/*
 * ************************************************************************
 *
 * Fill the component selector argument with the components needed to
 * allocate destination data.
 *
 * ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineSchedule::initializeDestinationVector(
   hier::ComponentSelector& allocate_vector) const
{
   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      allocate_vector.setFlag(d_refine_items[iri]->d_dst);
   }
}

}
}
