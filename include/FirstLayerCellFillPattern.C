//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Fill pattern class to provide interface for stencils
//

#ifndef included_pdat_FirstLayerCellFillPattern_C
#define included_pdat_FirstLayerCellFillPattern_C

#include "FirstLayerCellFillPattern.h"

#include "CellGeometry.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace pdat {

template<int DIM>
std::string FirstLayerCellFillPattern<DIM>::s_name_id =
   "FIRST_LAYER_CELL_FILL_PATTERN";

/*
*************************************************************************
*                                                                       *
* Constructor                                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
FirstLayerCellFillPattern<DIM>::FirstLayerCellFillPattern()
: d_stencil_width(1)
{
}

/*
*************************************************************************
*									*
* Destructor								*
*                                                                       *
*************************************************************************
*/
template<int DIM>
FirstLayerCellFillPattern<DIM>::~FirstLayerCellFillPattern()
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
FirstLayerCellFillPattern<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& dst_patch_box,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset) const
{
   NULL_USE(overwrite_interior);

   hier::BoxList<DIM> stencil_boxes;
   hier::Box<DIM> ghost_box(hier::Box<DIM>::grow(dst_patch_box,hier::IntVector<DIM>(1)));
   stencil_boxes.removeIntersections(ghost_box,dst_patch_box);

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

   hier::BoxOverlap<DIM> *overlap = new CellOverlap<DIM>(dst_boxes, src_offset);
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
hier::IntVector<DIM>& FirstLayerCellFillPattern<DIM>::getStencilWidth()
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
const std::string& FirstLayerCellFillPattern<DIM>::getPatternName() const
{
   return (s_name_id);
}

}
}
#endif
