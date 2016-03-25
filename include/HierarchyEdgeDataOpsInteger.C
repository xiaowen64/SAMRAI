//
// File:	HierarchyEdgeDataOpsInteger.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Operations for integer edge data on multiple levels.
//

#ifndef included_math_HierarchyEdgeDataOpsInteger_C
#define included_math_HierarchyEdgeDataOpsInteger_C

#include "HierarchyEdgeDataOpsInteger.h"
#include "PatchDescriptor.h"
#include "BoxUtilities.h"
#include "EdgeDataFactory.h"
#include "tbox/MPI.h"
#include "tbox/IEEE.h"
#include "tbox/Utilities.h"

#include <typeinfo>
using namespace std;
#include <stdlib.h>
#include <float.h>
#include <math.h>
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  HierarchyEdgeDataOpsInteger<DIM>::HierarchyEdgeDataOpsInteger(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
: HierarchyDataOpsInteger<DIM>()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
   if ( (coarsest_level < 0) || (finest_level < 0) ) {
      if ( d_hierarchy->getNumberLevels() == 0 ) {
         d_coarsest_level = coarsest_level;
         d_finest_level = finest_level;
      } else {
         resetLevels(0, d_hierarchy->getFinestLevelNumber());
      }
   } else {
      resetLevels(coarsest_level, finest_level);
   }
}
 
template<int DIM>  HierarchyEdgeDataOpsInteger<DIM>::~HierarchyEdgeDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Rotuines to set the hierarchy and level information.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::setPatchHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (coarsest_level >= 0)
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
            edge_boxes.getBox(i) =
               pdat::EdgeGeometry<DIM>::toEdgeBox(edge_boxes.getBox(i), nd);
         }
         hier::BoxUtilities<DIM>::makeNonOverlappingBoxLists(
                             d_nonoverlapping_edge_boxes[nd][ln],
                             edge_boxes);
      }
   }
}

template<int DIM> const tbox::Pointer< hier::PatchHierarchy<DIM> >
HierarchyEdgeDataOpsInteger<DIM>::getPatchHierarchy() const
{
   return(d_hierarchy);
}

/*
*************************************************************************
*                                                                       *
* Basic generic operations.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> int HierarchyEdgeDataOpsInteger<DIM>::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int entries = 0;

   if (interior_only) {

      tbox::Pointer< pdat::EdgeDataFactory<DIM,int> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < DIM; nd++) {
            assert(npatches == d_nonoverlapping_edge_boxes[nd][ln].getSize());
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
            tbox::Pointer< pdat::EdgeData<DIM,int> > d = 
               level->getPatch(ip())->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!d.isNull());
#endif
            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      int global_entries = tbox::MPI::sumReduction(entries);
      entries = global_entries;

   }

   return( entries );
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::copyData(
   const int dst_id, 
   const int src_id, 
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull()); 
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.copyData(d, s, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer< pdat::EdgeDataFactory<DIM,int> >
      d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   assert(!d1fact.isNull());
   tbox::Pointer< pdat::EdgeDataFactory<DIM,int> >
      d2fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id);
   assert(!d2fact.isNull());
   assert(d1fact->getDefaultDepth() == d2fact->getDefaultDepth());
   assert(d1fact->getDefaultGhostCellWidth() == d2fact->getDefaultGhostCellWidth());
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
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

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::printData(
   const int data_id,
   ostream& s,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   s << "Patch descriptor id = " << data_id << endl;
   s << "Factory = " << typeid(d_hierarchy->getPatchDescriptor()->
			      getPatchDataFactory(data_id)).name() << endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << endl;
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::setToScalar(
   const int data_id,
   const int& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
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

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::scale(
   const int dst_id,
   const int& alpha,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::addScalar(
   const int dst_id,
   const int src_id,
   const int& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::linearSum(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int& beta,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::axpy(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );
         
         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::axmy(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::EdgeData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.abs(d, src, box);
      }
   }
}

template<int DIM> int HierarchyEdgeDataOpsInteger<DIM>::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int minval = tbox::IEEE::getINT_MAX();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         minval = tbox::Utilities::imin( minval, d_patch_ops.min(d, box) );
      }
   }

   int global_min = tbox::MPI::minReduction(minval);
   return( global_min );
}

template<int DIM> int HierarchyEdgeDataOpsInteger<DIM>::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int maxval = -tbox::IEEE::getINT_MAX();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         maxval = tbox::Utilities::imax( maxval, d_patch_ops.min(d, box) );
      }
   }

   int global_max = tbox::MPI::maxReduction(maxval);
   return( global_max );
}

template<int DIM> void HierarchyEdgeDataOpsInteger<DIM>::setRandomValues(
   const int data_id,
   const int& width,
   const int& low,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::EdgeData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.setRandomValues(d, width, low, box);
      }
   }
}

}
}
#endif
