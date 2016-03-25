//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/standard/FillBoxSet.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Utility class for "smart" boxlist operations in comm schedules.
//

#ifndef included_xfer_FillBoxSet_C
#define included_xfer_FillBoxSet_C

#include "FillBoxSet.h"


#ifdef DEBUG_NO_INLINE
#include "FillBoxSet.I"
#endif

namespace SAMRAI {
   namespace xfer {

template<int DIM> 
FillBoxSet<DIM>::FillBoxSet(const hier::Box<DIM>& box)
{
   d_bounding_box = box;
   hier::BoxList<DIM> tmp(box);
   d_boxes = tmp;
   d_recompute_bounding_box = false;
}

template<int DIM> 
FillBoxSet<DIM>::FillBoxSet(
   const FillBoxSet<DIM>& fill_box_set) 
{
   d_bounding_box = fill_box_set.d_bounding_box;
   d_boxes = fill_box_set.d_boxes;
   d_recompute_bounding_box = false;
}

template<int DIM> 
FillBoxSet<DIM>::~FillBoxSet()
{
}

template<int DIM> 
FillBoxSet<DIM>& FillBoxSet<DIM>::operator=(
   const FillBoxSet<DIM>& fill_box_set) 
{
   if (this != &fill_box_set) {
      d_bounding_box = fill_box_set.d_bounding_box;
      d_boxes = fill_box_set.d_boxes;
      d_recompute_bounding_box = fill_box_set.d_recompute_bounding_box;
   }

   return(*this);
}

template<int DIM> 
void FillBoxSet<DIM>::resetFillBoxes(const hier::Box<DIM>& box)
{
   d_bounding_box = box;
   hier::BoxList<DIM> tmp(box);
   d_boxes = tmp;
   d_recompute_bounding_box = false;
}

template<int DIM> 
void FillBoxSet<DIM>::resetFillBoxes(const hier::BoxList<DIM>& boxes)
{
   d_bounding_box = boxes.getBoundingBox();
   d_boxes = boxes;
   d_recompute_bounding_box = false;
}

template<int DIM>
void FillBoxSet<DIM>::addFillBox(const hier::Box<DIM>& box)
{
   d_boxes.addItem(box);
   d_recompute_bounding_box = true;
}

template<int DIM> 
void FillBoxSet<DIM>::print(std::ostream& os) const
{
   os << "d_bounding_box = " << d_bounding_box << std::endl;
   os << "d_recompute_bounding_box = " << d_recompute_bounding_box << std::endl;
   os << "d_boxes = " << std::endl;
   d_boxes.print(os);
}

}
}

#endif
