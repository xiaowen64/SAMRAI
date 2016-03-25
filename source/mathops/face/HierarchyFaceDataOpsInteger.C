//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/face/HierarchyFaceDataOpsInteger.C $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Operations for integer face data on multiple levels.
//

#ifndef included_math_HierarchyFaceDataOpsInteger_C
#define included_math_HierarchyFaceDataOpsInteger_C

#include "HierarchyFaceDataOpsInteger.h"
#include "PatchDescriptor.h"
#include "BoxUtilities.h"
#include "FaceDataFactory.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"

#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
    namespace math {

template<int DIM>  HierarchyFaceDataOpsInteger<DIM>::HierarchyFaceDataOpsInteger(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
: HierarchyDataOpsInteger<DIM>()
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
 
template<int DIM>  HierarchyFaceDataOpsInteger<DIM>::~HierarchyFaceDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Rotuines to set the hierarchy and level information.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::setPatchHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::resetLevels(
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
   d_nonoverlapping_face_boxes[d].resizeArray(d_finest_level+1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      hier::BoxArray<DIM> face_boxes;
      const int n = level->getNumberOfPatches();

      for (int nd = 0; nd < DIM; nd++) {
         face_boxes = level->getBoxes();
         for (i = 0; i < n; i++) {
            face_boxes[i] =
               pdat::FaceGeometry<DIM>::toFaceBox(face_boxes[i], nd);
         }
         hier::BoxUtilities<DIM>::makeNonOverlappingBoxLists(
                             d_nonoverlapping_face_boxes[nd][ln],
                             face_boxes);
      }
   }
}

template<int DIM> const tbox::Pointer< hier::PatchHierarchy<DIM> >
HierarchyFaceDataOpsInteger<DIM>::getPatchHierarchy() const
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

template<int DIM> int HierarchyFaceDataOpsInteger<DIM>::numberOfEntries(
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

      tbox::Pointer< pdat::FaceDataFactory<DIM,int> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < DIM; nd++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_face_boxes[nd][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {
            typename tbox::List< hier::Box<DIM> >::Iterator lb;
            for (int eb = 0; eb < DIM; eb++) {
               lb = ((d_nonoverlapping_face_boxes[eb][ln])[il]).listStart();
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
            tbox::Pointer< pdat::FaceData<DIM,int> > d = 
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

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::copyData(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull()); 
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.copyData(d, s, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer< pdat::FaceDataFactory<DIM,int> >
      d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   TBOX_ASSERT(!d1fact.isNull());
   tbox::Pointer< pdat::FaceDataFactory<DIM,int> >
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

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::printData(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::setToScalar(
   const int data_id,
   const int& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(data_id);
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

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::scale(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::addScalar(
   const int dst_id,
   const int src_id,
   const int& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::add(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::subtract(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::multiply(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::divide(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::reciprocal(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::linearSum(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int& beta,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::axpy(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );
         
         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::axmy(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::abs(
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.abs(d, src, box);
      }
   }
}

template<int DIM> int HierarchyFaceDataOpsInteger<DIM>::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int minval = tbox::MathUtilities<int>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         minval = tbox::MathUtilities<int>::Min( minval, 
                                                 d_patch_ops.min(d, box) );
      }
   }

   int global_min = tbox::SAMRAI_MPI::minReduction(minval);
   return( global_min );
}

template<int DIM> int HierarchyFaceDataOpsInteger<DIM>::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int maxval = -(tbox::MathUtilities<int>::getMax());

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         maxval = tbox::MathUtilities<int>::Max( maxval, 
                                                 d_patch_ops.min(d, box) );
      }
   }

   int global_max = tbox::SAMRAI_MPI::maxReduction(maxval);
   return( global_max );
}

template<int DIM> void HierarchyFaceDataOpsInteger<DIM>::setRandomValues(
   const int data_id,
   const int& width,
   const int& low,
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

         tbox::Pointer< pdat::FaceData<DIM,int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.setRandomValues(d, width, low, box);
      }
   }
}

}
}
#endif
