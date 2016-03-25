//
// File:	SkeletonRefine.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Constant refine operator for cell-centered double data on 
//              a Moving mesh.
//

#ifndef included_geom_SkeletonRefine_C
#define included_geom_SkeletonRefine_C

#include "SkeletonRefine.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace geom {

template<int DIM> SkeletonRefine<DIM>::SkeletonRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "SKELETON_REFINE";
}

template<int DIM> SkeletonRefine<DIM>::~SkeletonRefine()
{
}

template<int DIM> bool SkeletonRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   NULL_USE(var);
   if (op_name == d_name_id) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
SkeletonRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int SkeletonRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
SkeletonRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void SkeletonRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   //no operation for the empty refine operator
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(dst_component);
   NULL_USE(src_component);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

}
}
#endif
