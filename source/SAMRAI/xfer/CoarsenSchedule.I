/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Coarsening schedule for data transfer between AMR levels
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 * ************************************************************************
 *
 * Return const pointer to equivalence classes used in schedule.
 *
 * ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<CoarsenClasses>&
CoarsenSchedule::getEquivalenceClasses() const
{
   return d_coarsen_classes;
}

}
}
