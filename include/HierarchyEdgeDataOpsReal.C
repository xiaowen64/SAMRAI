//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/edge/HierarchyEdgeDataOpsReal.C $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Templated operations for real edge data on multiple levels.
//

#ifndef included_math_HierarchyEdgeDataOpsReal_C
#define included_math_HierarchyEdgeDataOpsReal_C

#include "HierarchyEdgeDataOpsReal.h"
#include "BoxUtilities.h"
#include "PatchDescriptor.h"
#include "EdgeDataFactory.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
HierarchyEdgeDataOpsReal<DIM,TYPE>::HierarchyEdgeDataOpsReal(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
: HierarchyDataOpsReal<DIM,TYPE>()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
   if ( (coarsest_level < 0) || (finest_level < 0) ) {
      if ( d_hierarchy->getNumberOfLevels() == 0 ) {
         d_coarsest_level = coarsest_level;
         d_finest_level = finest_level;
      } else {
         resetLevels(0, d_hierarchy->getFinestLevelNumber());
      }
   } else {
      resetLevels(coarsest_level, finest_level);
   }
}
 
template<int DIM, class TYPE>
HierarchyEdgeDataOpsReal<DIM,TYPE>::~HierarchyEdgeDataOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* Routines to set the hierarchy and level informtation.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::setPatchHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (coarsest_level >= 0)
          && (finest_level >= coarsest_level)
          && (finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   for (int d = 0; d < DIM; d++) {
      d_nonoverlapping_edge_boxes[d].resizeArray(d_finest_level+1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      hier::BoxArray<DIM> edge_boxes;
      const int n = level->getNumberOfPatches();

      for (int nd = 0; nd < DIM; nd++ ) {
         edge_boxes = level->getBoxes();
         for (i = 0; i < n; i++) {
            edge_boxes[i] =
               pdat::EdgeGeometry<DIM>::toEdgeBox(edge_boxes[i], nd);
         }
         hier::BoxUtilities<DIM>::makeNonOverlappingBoxLists(
                             d_nonoverlapping_edge_boxes[nd][ln],
                             edge_boxes);
      }
   }
}

template<int DIM, class TYPE>
const tbox::Pointer< hier::PatchHierarchy<DIM> >
HierarchyEdgeDataOpsReal<DIM,TYPE>::getPatchHierarchy() const
{
   return(d_hierarchy);
}

/*
*************************************************************************
*                                                                       *
* The following are private and cannot be used, but they are defined    *
* here for compilers that require that every template declaration have  *
* a definition (a stupid requirement, if you ask me).                   *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
HierarchyEdgeDataOpsReal<DIM,TYPE>::HierarchyEdgeDataOpsReal(
   const HierarchyEdgeDataOpsReal<DIM,TYPE>& foo)
:  HierarchyDataOpsReal<DIM,TYPE>()
{
   NULL_USE(foo); 
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::operator=(
   const HierarchyEdgeDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

/*
*************************************************************************
*                                                                       *
* Basic generic operations.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::copyData(
   const int dst_id, 
   const int src_id, 
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull()); 
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.copyData(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer< pdat::EdgeDataFactory<DIM,TYPE> >
      d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   TBOX_ASSERT(!d1fact.isNull());
   tbox::Pointer< pdat::EdgeDataFactory<DIM,TYPE> >
      d2fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id);
   TBOX_ASSERT(!d2fact.isNull());
   TBOX_ASSERT(d1fact->getDefaultDepth() == d2fact->getDefaultDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   s << "Patch descriptor id = " << data_id << std::endl;
   s << "Factory = " << typeid(*d_hierarchy->getPatchDescriptor()->
			       getPatchDataFactory(data_id)).name() << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << std::endl;
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::setToScalar(
   const int data_id,
   const TYPE& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.setToScalar(d, alpha, box);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Basic generic arithmetic operations.                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::scale(
   const int dst_id,
   const TYPE& alpha,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::addScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.add(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.subtract(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.multiply(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.divide(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.reciprocal(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::linearSum(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const TYPE& beta,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.linearSum(dst, alpha, src1, beta, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::axpy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );
         
         d_patch_ops.axpy(dst, alpha, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::axmy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.axmy(dst, alpha, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.abs(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::setRandomValues(
   const int data_id,
   const TYPE& width,
   const TYPE& low,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!data.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : data->getGhostBox() );

         d_patch_ops.setRandomValues(data, width, low, box);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Generic norm and order operations.                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int HierarchyEdgeDataOpsReal<DIM,TYPE>::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int entries = 0;

   if (interior_only) {

      tbox::Pointer< pdat::EdgeDataFactory<DIM,TYPE> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < DIM; nd++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_edge_boxes[nd][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {
            typename tbox::List< hier::Box<DIM> >::Iterator lb;

            for (int eb = 0; eb < DIM; eb++) {
               lb = ((d_nonoverlapping_edge_boxes[eb][ln])[il]).listStart();
               for ( ;lb;lb++) {
                  entries += lb().size();
               }
            }
         }
      }

      entries *= dfact->getDefaultDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
            tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d = 
               level->getPatch(ip())->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      int global_entries = tbox::SAMRAI_MPI::sumReduction(entries);
      entries = global_entries;

   }

   return( entries );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(vol_id >= 0);
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data  = p->getPatchData(data_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv  = p->getPatchData(vol_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!cv.isNull());
#endif
         hier::Box<DIM> box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(data, cv, box);
      }
   }

   double global_sum = tbox::SAMRAI_MPI::sumReduction(sum);
   return( global_sum );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::L1Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm += d_patch_ops.L1Norm(data, box, cv);
      }
   }

   if ( !local_only ) {
      norm = tbox::SAMRAI_MPI::sumReduction(norm);
   }
   return( norm );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::L2Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   double norm_squared = HierarchyEdgeDataOpsReal<DIM,TYPE>::dot(data_id,
                                                                   data_id,
                                                                   vol_id,
                                                                   local_only);

   return( sqrt(norm_squared) );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::weightedL2Norm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm_squared = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data   = p->getPatchData(data_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > weight = p->getPatchData(wgt_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         double pnorm = d_patch_ops.weightedL2Norm(data, weight, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = tbox::SAMRAI_MPI::sumReduction(norm_squared);
   return( sqrt(global_norm_squared) );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);
   
   double volume = ( (vol_id < 0) ? (double) numberOfEntries(data_id, true)
                                  : sumControlVolumes(data_id, vol_id) );

   double rms_norm = l2_norm/sqrt(volume);
   return( rms_norm );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::weightedRMSNorm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{

   double l2_norm = weightedL2Norm(data_id, wgt_id, vol_id);

   double volume = ( (vol_id < 0) ? (double) numberOfEntries(data_id, true)
                                  : sumControlVolumes(data_id, vol_id) );

   double rms_norm = l2_norm/sqrt(volume);
   return( rms_norm );
}

template<int DIM, class TYPE>
double HierarchyEdgeDataOpsReal<DIM,TYPE>::maxNorm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm = tbox::MathUtilities<double>::Max(norm, 
                                            d_patch_ops.maxNorm(data, box, cv));
      }
   }

   if ( !local_only ) {
      norm = tbox::SAMRAI_MPI::maxReduction(norm);
   }
   return( norm );
}

template<int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE dprod = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data1.isNull());
#endif
            box = data1->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         dprod += d_patch_ops.dot(data1, data2, box, cv);
      }
   }

   if ( local_only ) {
      return dprod;
   }
   else {
      TYPE global_dot = tbox::SAMRAI_MPI::sumReduction(dprod);
      return( global_dot );
   }
}


template<int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::integral(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE local_integral = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data =
                                                  p->getPatchData(data_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > vol = p->getPatchData(vol_id);

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!data.isNull());
         TBOX_ASSERT(!vol.isNull());
#endif

         hier::Box<DIM> box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   TYPE global_integral = tbox::SAMRAI_MPI::sumReduction(local_integral);
   return( global_integral );
}


/*
*************************************************************************
*                                                                       *
* Generic miscellaneous operations for real data.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int HierarchyEdgeDataOpsReal<DIM,TYPE>::computeConstrProdPos(
   const int data1_id, 
   const int data2_id, 
   const int vol_id) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data1.isNull());
#endif
            box = data1->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         test = tbox::MathUtilities<int>::Min( test, 
                   d_patch_ops.computeConstrProdPos(data1, data2, box, cv) );
      }
   }

   int global_test = tbox::SAMRAI_MPI::minReduction(test);
   return( global_test );
}

template<int DIM, class TYPE>
void HierarchyEdgeDataOpsReal<DIM,TYPE>::compareToScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src  = p->getPatchData(src_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }
  
         d_patch_ops.compareToScalar(dst, src, alpha, box, cv);
      }
   }
}

template<int DIM, class TYPE>
int HierarchyEdgeDataOpsReal<DIM,TYPE>::testReciprocal(
   const int dst_id, 
   const int src_id, 
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > src  = p->getPatchData(src_id);
         tbox::Pointer< pdat::EdgeData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         test = tbox::MathUtilities<int>::Min( test,
                   d_patch_ops.testReciprocal(dst, src, box, cv) );
      }
   }

   int global_test = tbox::SAMRAI_MPI::minReduction(test);
   return( global_test );
}

template<int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const int numer_id, 
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE max = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > numer  = p->getPatchData(numer_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > denom  = p->getPatchData(denom_id);

         hier::Box<DIM> box = p->getBox();

         max = tbox::MathUtilities<TYPE>::Max(max,
                   d_patch_ops.maxPointwiseDivide(numer, denom, box));
      }
   }

   if ( !local_only ) {
      max = tbox::SAMRAI_MPI::maxReduction(max);
   }
   return( max );
}


template <int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::minPointwiseDivide(
   const int numer_id, 
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE min = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > numer  = p->getPatchData(numer_id);
         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > denom  = p->getPatchData(denom_id);

         hier::Box<DIM> box = p->getBox();

         min = tbox::MathUtilities<TYPE>::Min(min,
                   d_patch_ops.minPointwiseDivide(numer, denom, box));
      }
   }

   if ( !local_only ) {
      min = tbox::SAMRAI_MPI::minReduction(min);
   }
   return( min );
}

template<int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         minval = tbox::MathUtilities<TYPE>::Min( minval, d_patch_ops.min(d, box) );
      }
   }

   TYPE global_min = tbox::SAMRAI_MPI::minReduction(minval);
   return( global_min );
}

template<int DIM, class TYPE>
TYPE HierarchyEdgeDataOpsReal<DIM,TYPE>::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         maxval = tbox::MathUtilities<TYPE>::Max( maxval, d_patch_ops.max(d, box) );
      }
   }

   TYPE global_max = tbox::SAMRAI_MPI::maxReduction(maxval);
   return( global_max );
}

}
}
#endif
