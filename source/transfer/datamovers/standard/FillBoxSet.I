//
// File:	FillBoxSet.I
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Utility class for "smart" boxlist operations in comm schedules.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace xfer {

template<int DIM> inline
FillBoxSet<DIM>::FillBoxSet()
{
   d_recompute_bounding_box = false;
}

template<int DIM> inline
const hier::Box<DIM>& FillBoxSet<DIM>::getBoundingBox()
{
   if (d_recompute_bounding_box) {
      d_bounding_box = d_boxes.getBoundingBox();
      d_recompute_bounding_box = false;
   }
   return(d_bounding_box);
}

template<int DIM> inline
const hier::BoxList<DIM>& FillBoxSet<DIM>::getBoxList() const
{
   return(d_boxes);
}

template<int DIM> inline
hier::BoxList<DIM>& FillBoxSet<DIM>::getBoxListToChange()
{
   d_recompute_bounding_box = true;
   return(d_boxes);
}

template<int DIM> inline
int FillBoxSet<DIM>::getNumberOfBoxes() const
{
   return(d_boxes.getNumberOfBoxes());
}

template<int DIM> inline
void FillBoxSet<DIM>::simplifyBoxes()
{
   d_boxes.simplifyBoxes();
}

template<int DIM> inline
void FillBoxSet<DIM>::removeIntersections(const hier::Box<DIM>& takeaway)
{
   if (getBoundingBox().intersects(takeaway)) {
      d_boxes.removeIntersections(takeaway);
      d_bounding_box = d_boxes.getBoundingBox();
      d_recompute_bounding_box = false;
   }
}

template<int DIM> inline
void FillBoxSet<DIM>::removeIntersections(const hier::BoxList<DIM>& takeaway)
{
   d_boxes.removeIntersections(takeaway);
   d_bounding_box = d_boxes.getBoundingBox();
   d_recompute_bounding_box = false;
}

template<int DIM> inline
void FillBoxSet<DIM>::intersectBoxes(const hier::Box<DIM>& box)
{
   if (d_bounding_box.intersects(box)) {
      d_boxes.intersectBoxes(box);
      d_bounding_box = d_boxes.getBoundingBox();
      d_recompute_bounding_box = false;
   }
}

template<int DIM> inline
void FillBoxSet<DIM>::intersectBoxes(const hier::BoxList<DIM>& boxes)
{
   d_boxes.intersectBoxes(boxes);
   d_bounding_box = d_boxes.getBoundingBox();
   d_recompute_bounding_box = false;
}

#ifdef DEBUG_NO_TEMPLATE 
#undef inline inline
#endif

}
}
