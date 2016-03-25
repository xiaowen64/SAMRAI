//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Abstract fill pattern class to provide interface for stencils
//

#ifndef included_xfer_BoxGeometryFillPattern_C
#define included_xfer_BoxGeometryFillPattern_C

#include "BoxGeometryFillPattern.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace xfer {


template<int DIM> std::string BoxGeometryFillPattern<DIM>::s_name_id = "BOX_GEOMETRY_FILL_PATTERN";
template<int DIM>
hier::IntVector<DIM> BoxGeometryFillPattern<DIM>::s_stencil_width =
   hier::IntVector<DIM>(0);

/*
*************************************************************************
*                                                                       *
* Default contructor only sets the string name identifier               *
*                                                                       *
*************************************************************************
*/
template<int DIM>
BoxGeometryFillPattern<DIM>::BoxGeometryFillPattern()
{
}

/*
*************************************************************************
*									*
* Destructor                                                            *
*       								*
*************************************************************************
*/
template<int DIM>
BoxGeometryFillPattern<DIM>::~BoxGeometryFillPattern()
{
}

/*
*************************************************************************
*                                                                       *
* Calculate the overlap using the implemented calculateOverlap() method *
* for the destination geometry.                                         *
*                                                                       *
*************************************************************************
*/
template<int DIM>
tbox::Pointer< hier::BoxOverlap<DIM> >
BoxGeometryFillPattern<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& dst_patch_box,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset) const
{
   NULL_USE(dst_patch_box);
   return (dst_geometry.calculateOverlap(src_geometry, src_mask,
                                         overwrite_interior, src_offset) );
}

/*
*************************************************************************
*                                                                       *
* Return the string name identifier.                                    *
*                                                                       *   
*************************************************************************
*/
template<int DIM>
const std::string& BoxGeometryFillPattern<DIM>::getPatternName() const
{
   return (s_name_id);
}

/*
*************************************************************************
*                                                                       *
* getStencilWidth() throws an error if called.  Only overridding        *
* versions of this method in concrete subclasses should be called.      *
*                                                                       *
*************************************************************************
*/
template<int DIM>
hier::IntVector<DIM>& BoxGeometryFillPattern<DIM>::getStencilWidth()
{
   TBOX_ERROR("BoxGeometryFillPattern<DIM>::getStencilWidth() should not be\n"
              << "called.  This pattern creates overlaps based on\n"
              << "the BoxGeometry objects and is not restricted to a\n"
              << "specific stencil.\n");

   return (s_stencil_width);
}


}
}
#endif
