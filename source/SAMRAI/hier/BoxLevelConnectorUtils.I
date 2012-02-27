/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utilities for working on DLBG edges.
 *
 ************************************************************************/
namespace SAMRAI {
namespace hier {

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::setSanityCheckMethodPreconditions(
   bool do_check)
{
   d_sanity_check_precond = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::setSanityCheckMethodPostconditions(
   bool do_check)
{
   d_sanity_check_postcond = do_check;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::computeExternalParts(
   BoxLevel& external,
   Connector& input_to_external,
   const Connector& input_to_reference,
   const IntVector& nesting_width,
   const BoxContainer& domain) const
{
   t_compute_external_parts->start();

   computeInternalOrExternalParts(
      external,
      input_to_external,
      'e',
      input_to_reference,
      nesting_width,
      domain);

   t_compute_external_parts->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::computeInternalParts(
   BoxLevel& internal,
   Connector& input_to_internal,
   const Connector& input_to_reference,
   const IntVector& nesting_width,
   const BoxContainer& domain) const
{
   t_compute_internal_parts->start();

   computeInternalOrExternalParts(
      internal,
      input_to_internal,
      'i',
      input_to_reference,
      nesting_width,
      domain);

   t_compute_internal_parts->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::initializeCallback()
{
   t_make_sorting_map = tbox::TimerManager::getManager()->
      getTimer("BoxLevelConnectorUtils::makeSortingMap()");
   t_compute_external_parts = tbox::TimerManager::getManager()->
      getTimer("BoxLevelConnectorUtils::computeExternalParts()");
   t_compute_external_parts_intersection =
      tbox::TimerManager::getManager()->
      getTimer("BoxLevelConnectorUtils::computeExternalParts()_intersection");
   t_compute_internal_parts = tbox::TimerManager::getManager()->
      getTimer("BoxLevelConnectorUtils::computeInternalParts()");
   t_compute_internal_parts_intersection =
      tbox::TimerManager::getManager()->
      getTimer("BoxLevelConnectorUtils::computeInternalParts()_intersection");
}

/*
 ***************************************************************************
 * * Release static timers.  To be called by shutdown registry to make sure
 * * memory for timers does not leak.
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
BoxLevelConnectorUtils::finalizeCallback()
{
   t_make_sorting_map.reset();
   t_compute_external_parts.reset();
   t_compute_external_parts_intersection.reset();
   t_compute_internal_parts.reset();
   t_compute_internal_parts_intersection.reset();
}

}
}
