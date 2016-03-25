//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/geometry/skeleton/operators/SkeletonCoarsen.C $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Weighted averaging operator for cell-centered double data on 
//              a Moving mesh.
//

#ifndef included_geom_SkeletonCoarsen_C
#define included_geom_SkeletonCoarsen_C

#include "SkeletonCoarsen.h"

#include "tbox/Utilities.h"

#include<float.h>
#include<math.h>


namespace SAMRAI {
    namespace geom {

template<int DIM> SkeletonCoarsen<DIM>::SkeletonCoarsen()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "SKELETON_COARSEN";
}

template<int DIM> SkeletonCoarsen<DIM>::~SkeletonCoarsen()
{
}

template<int DIM> bool SkeletonCoarsen<DIM>::findCoarsenOperator(
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

template<int DIM> const std::string& SkeletonCoarsen<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int SkeletonCoarsen<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> SkeletonCoarsen<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void SkeletonCoarsen<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   //no operation for skeleton coarsen operator
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(dst_component);
   NULL_USE(src_component);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

}
}
#endif
