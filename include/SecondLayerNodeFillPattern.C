//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Fill pattern class to provide interface for stencils
//

#ifndef included_pdat_SecondLayerNodeFillPattern_C
#define included_pdat_SecondLayerNodeFillPattern_C

#include "SecondLayerNodeFillPattern.h"

#include "NodeGeometry.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace pdat {

template<int DIM>
std::string SecondLayerNodeFillPattern<DIM>::s_name_id =
   "SECOND_LAYER_NODE_FILL_PATTERN";

/*
*************************************************************************
*                                                                       *
* Constructor                                                           *
*                                                                       *
*************************************************************************
*/
template<int DIM>
SecondLayerNodeFillPattern<DIM>::SecondLayerNodeFillPattern()
: d_stencil_width(0)
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
SecondLayerNodeFillPattern<DIM>::~SecondLayerNodeFillPattern()
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
SecondLayerNodeFillPattern<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& dst_patch_box,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset) const
{
   NULL_USE(overwrite_interior);

   hier::Box<DIM> dst_node_box(pdat::NodeGeometry<DIM>::toNodeBox(dst_patch_box));
   hier::Box<DIM> src_node_mask(pdat::NodeGeometry<DIM>::toNodeBox(src_mask));
   bool corner_overlap = ((dst_node_box * src_node_mask).size() == 1)
                         ? true : false;

   hier::BoxList<DIM> stencil_boxes;
   hier::Box<DIM> ghost_box(dst_node_box);
   ghost_box.grow(hier::IntVector<DIM>(1));
   stencil_boxes.removeIntersections(ghost_box, dst_node_box);
   if (corner_overlap) {
      hier::IntVector<DIM> grow_vec(0);
      hier::Box<DIM> grow_box;
      hier::BoxList<DIM> remove_list;
      for (unsigned int i = 0; i < DIM; i++) {
         grow_box = dst_node_box;
         grow_vec(i) = 1;
         grow_box.grow(grow_vec);
         remove_list.addItem(grow_box);
         grow_vec(i) = 0;
      }
      stencil_boxes.removeIntersections(remove_list);
   } 

   hier::BoxList<DIM> dst_boxes;

   const NodeGeometry<DIM> *t_dst =
      dynamic_cast<const NodeGeometry<DIM> *>(&dst_geometry);
   const NodeGeometry<DIM> *t_src =
      dynamic_cast<const NodeGeometry<DIM> *>(&src_geometry);

   TBOX_ASSERT(t_dst);
   TBOX_ASSERT(t_src);

   t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask,
                                  false, src_offset);

   dst_boxes.intersectBoxes(stencil_boxes);

   hier::BoxOverlap<DIM> *overlap = new NodeOverlap<DIM>(dst_boxes, src_offset);
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
hier::IntVector<DIM>& SecondLayerNodeFillPattern<DIM>::getStencilWidth()
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
const std::string& SecondLayerNodeFillPattern<DIM>::getPatternName() const
{
   return (s_name_id);
}



}
}
#endif
