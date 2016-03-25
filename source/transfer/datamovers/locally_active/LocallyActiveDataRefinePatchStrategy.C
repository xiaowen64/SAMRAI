//
// File:	LocallyActiveDataRefinePatchStrategy.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 306 $
// Modified:	$Date: 2005-04-26 13:49:29 -0700 (Tue, 26 Apr 2005) $
// Description:	Strategy interface to user routines for refining locally-active AMR data.
//

#ifndef included_xfer_LocallyActiveRefinePatchStrategy_C
#define included_xfer_LocallyActiveRefinePatchStrategy_C

#include "LocallyActiveDataRefinePatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*									*
* The default constructor and virtual destructor do nothing             *
* particularly interesting.		                                *
*									*
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefinePatchStrategy<DIM>::LocallyActiveDataRefinePatchStrategy()
{
}

template<int DIM>
LocallyActiveDataRefinePatchStrategy<DIM>::~LocallyActiveDataRefinePatchStrategy()
{
}

/*
*************************************************************************
*									*
* Loop over all fill boxes and call the user-defined preprocesses.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefinePatchStrategy<DIM>::setPhysicalBoundaryConditions(
   hier::Patch<DIM>& patch,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& refine_data,
   const double fill_time,
   const hier::IntVector<DIM>& ghost_width_to_fill)
{
   tbox::List<int> scratch_data_ids;
   for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
      vi(refine_data); vi; vi++) {
      scratch_data_ids.appendItem(vi()->d_scratch);
   }
   this->setPhysicalBoundaryConditions(patch, scratch_data_ids,
                                       fill_time, ghost_width_to_fill);
}

/*
*************************************************************************
*									*
* Loop over all fill boxes and call the user-defined preprocesses.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefinePatchStrategy<DIM>::preprocessRefineBoxes(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const xfer::LocallyActiveDataFillBoxSet<DIM>& fine_boxes,
   const hier::IntVector<DIM>& ratio)
{
   typename tbox::List< typename xfer::LocallyActiveDataFillBox<DIM> >::Iterator
      fbi(fine_boxes.getLocallyActiveDataBoxes());
   for ( ; fbi; fbi++) {
      tbox::List<int> scratch_data_ids;
      for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
           vi(fbi().getActiveRefineVarData()); vi; vi++) {
         scratch_data_ids.appendItem(vi()->d_scratch);
      }
      this->preprocessRefine(fine, coarse, scratch_data_ids,
                             fbi().getBox(), ratio);
      scratch_data_ids.clearItems();
   }
}

/*
*************************************************************************
*									*
* Loop over all fill boxes and call the user-defined postprocesses.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefinePatchStrategy<DIM>::postprocessRefineBoxes(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const xfer::LocallyActiveDataFillBoxSet<DIM>& fine_boxes,
   const hier::IntVector<DIM>& ratio)
{
   typename tbox::List< typename xfer::LocallyActiveDataFillBox<DIM> >::Iterator
      fbi(fine_boxes.getLocallyActiveDataBoxes());
   for ( ; fbi; fbi++) {
      tbox::List<int> scratch_data_ids;
      for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
           vi(fbi().getActiveRefineVarData()); vi; vi++) {
         scratch_data_ids.appendItem(vi()->d_scratch);
      }
      this->postprocessRefine(fine, coarse, scratch_data_ids,
                              fbi().getBox(), ratio);
      scratch_data_ids.clearItems();
   }
}

}
}
#endif
