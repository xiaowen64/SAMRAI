//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Fill pattern class to provide interface for stencils
//

#ifndef included_pdat_FirstLayerCellNoCornersFillPattern_C
#define included_pdat_FirstLayerCellNoCornersFillPattern_C

#include "FirstLayerCellNoCornersFillPattern.h"

#include "CellGeometry.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace pdat {

template<int DIM> std::string FirstLayerCellNoCornersFillPattern<DIM>::s_name_id =
   "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN";

/*
*************************************************************************
*                                                                       *
* Constructor                                                           *
*                                                                       *
*************************************************************************
*/
template<int DIM>
FirstLayerCellNoCornersFillPattern<DIM>::FirstLayerCellNoCornersFillPattern()
: d_stencil_width(1)
{
}

/*
*************************************************************************
*                                                                       *
* Destructor                                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM>
FirstLayerCellNoCornersFillPattern<DIM>::~FirstLayerCellNoCornersFillPattern()
{
}

/*
*************************************************************************
*                                                                       *
* Calculate the overlap according to the desired pattern                *
*                                                                       *
*************************************************************************
*/
template<int DIM>
tbox::Pointer< hier::BoxOverlap<DIM> >
FirstLayerCellNoCornersFillPattern<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& dst_patch_box,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset) const
{
   NULL_USE(overwrite_interior);

   hier::BoxList<DIM> stencil_boxes;
   for (unsigned short i = 0; i < DIM; i++) {
      hier::Box<DIM> low_box(dst_patch_box);
      low_box.lower(i) = dst_patch_box.lower(i) - 1;
      low_box.upper(i) = low_box.lower(i);
      stencil_boxes.addItem(low_box); 

      hier::Box<DIM> high_box(dst_patch_box);
      high_box.lower(i) = dst_patch_box.upper(i) + 1;
      high_box.upper(i) = high_box.lower(i);
      stencil_boxes.addItem(high_box); 
   }

   hier::BoxList<DIM> dst_boxes;

   const CellGeometry<DIM> *t_dst =
      dynamic_cast<const CellGeometry<DIM> *>(&dst_geometry);
   const CellGeometry<DIM> *t_src =
      dynamic_cast<const CellGeometry<DIM> *>(&src_geometry);

   TBOX_ASSERT(t_dst);
   TBOX_ASSERT(t_src);
   t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask,
                                  false, src_offset);

   dst_boxes.intersectBoxes(stencil_boxes);

   hier::BoxOverlap<DIM> *overlap =
      new CellOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));

}

/*
*************************************************************************
*                                                                       *
* Return the stencil width (1)                                          *
*                                                                       *
*************************************************************************
*/
template<int DIM>
hier::IntVector<DIM>& FirstLayerCellNoCornersFillPattern<DIM>::getStencilWidth()
{
   return (d_stencil_width);
}

/*
*************************************************************************
*                                                                       *
* Return the string name identifier                                     *
*                                                                       *
*************************************************************************
*/
template<int DIM>
const std::string& FirstLayerCellNoCornersFillPattern<DIM>::getPatternName() const
{
   return (s_name_id);
}


}
}
#endif
