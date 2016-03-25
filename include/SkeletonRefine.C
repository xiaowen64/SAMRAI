//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/skeleton/operators/SkeletonRefine.C $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
   const std::string &op_name) const
{
   NULL_USE(var);
   if (op_name == d_name_id) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
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
