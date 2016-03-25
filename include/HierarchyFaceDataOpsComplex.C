//
// File:	HierarchyFaceDataOpsComplex.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Operations for complex face data on multiple levels.
//

#ifndef included_math_HierarchyFaceDataOpsComplex_C
#define included_math_HierarchyFaceDataOpsComplex_C

#include "HierarchyFaceDataOpsComplex.h"
#include "BoxUtilities.h"
#include "PatchDescriptor.h"
#include "FaceDataFactory.h"
#include "tbox/MPI.h"
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

template<int DIM>  HierarchyFaceDataOpsComplex<DIM>::HierarchyFaceDataOpsComplex(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
: HierarchyDataOpsComplex<DIM>()
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
 
template<int DIM>  HierarchyFaceDataOpsComplex<DIM>::~HierarchyFaceDataOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Routines to set the hierarchy and level information.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::setPatchHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::resetLevels(
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
   d_nonoverlapping_face_boxes[d].resizeArray(d_finest_level+1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      hier::BoxArray<DIM> face_boxes;
      const int n = level->getNumberOfPatches();

      for (int nd = 0; nd < DIM; nd++) {
         face_boxes = level->getBoxes();
         for (i = 0; i < n; i++) {
            face_boxes.getBox(i) =
               pdat::FaceGeometry<DIM>::toFaceBox(face_boxes.getBox(i), nd);
         }
         hier::BoxUtilities<DIM>::makeNonOverlappingBoxLists(
                             d_nonoverlapping_face_boxes[nd][ln],
                             face_boxes);
      }
   }
}

template<int DIM> const tbox::Pointer< hier::PatchHierarchy<DIM> >
HierarchyFaceDataOpsComplex<DIM>::getPatchHierarchy() const
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

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::copyData(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull()); 
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.copyData(d, s, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer< pdat::FaceDataFactory<DIM,dcomplex> >
      d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   assert(!d1fact.isNull());
   tbox::Pointer< pdat::FaceDataFactory<DIM,dcomplex> >
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

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::printData(
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
   s << "Factory = " << typeid(*d_hierarchy->getPatchDescriptor()->
			       getPatchDataFactory(data_id)).name() << endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << endl;
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::setToScalar(
   const int data_id,
   const dcomplex& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(data_id);
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

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::scale(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::addScalar(
   const int dst_id,
   const int src_id,
   const dcomplex& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::add(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::subtract(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::multiply(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::divide(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::reciprocal(
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::linearSum(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const dcomplex& beta,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::axpy(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );
         
         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::axmy(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::abs(
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

         tbox::Pointer< pdat::FaceData<DIM,double> > d = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.abs(d, src, box);
      }
   }
}

template<int DIM> void HierarchyFaceDataOpsComplex<DIM>::setRandomValues(
   const int data_id,
   const dcomplex& width,
   const dcomplex& low,
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

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.setRandomValues(d, width, low, box);
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

template<int DIM> int HierarchyFaceDataOpsComplex<DIM>::numberOfEntries(
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

      tbox::Pointer< pdat::FaceDataFactory<DIM,dcomplex> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < DIM; nd++) {
            assert(npatches == d_nonoverlapping_face_boxes[nd][ln].getSize());
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
            tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = 
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

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(vol_id >= 0);
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv  = p->getPatchData(vol_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!cv.isNull());
#endif
         hier::Box<DIM> box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(d, cv, box);
      }
   }

   double global_sum = tbox::MPI::sumReduction(sum);
   return( global_sum );
}

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::L1Norm(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d  = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!d.isNull());
#endif
            box = d->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm += d_patch_ops.L1Norm(d, box, cv);
      }
   }

   double global_norm = tbox::MPI::sumReduction(norm);
   return( global_norm );
}

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::L2Norm(
   const int data_id,
   const int vol_id) const
{
   dcomplex dotprod = HierarchyFaceDataOpsComplex<DIM>::dot(data_id,
                                                             data_id,
                                                             vol_id);

   return( sqrt(real(dotprod)) );
}

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::weightedL2Norm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm_squared = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > w = p->getPatchData(wgt_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!d.isNull());
#endif
            box = d->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         double pnorm = d_patch_ops.weightedL2Norm(d, w, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = tbox::MPI::sumReduction(norm_squared);
   return( sqrt(global_norm_squared) );
}

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);
   
   double volume = ( (vol_id < 0) ? (double) numberOfEntries(data_id, true)
                                  : sumControlVolumes(data_id, vol_id) );

   double rms_norm = l2_norm/sqrt(volume);
   return( rms_norm );
}

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::weightedRMSNorm(
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

template<int DIM> double HierarchyFaceDataOpsComplex<DIM>::maxNorm(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!d.isNull());
#endif
            box = d->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm = tbox::Utilities::dmax(norm, d_patch_ops.maxNorm(d, box, cv));
      }
   }

   double global_norm = tbox::MPI::maxReduction(norm);
   return( global_norm );
}

template<int DIM> dcomplex HierarchyFaceDataOpsComplex<DIM>::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   dcomplex dprod = dcomplex(0.0,0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d1 = 
                                                  p->getPatchData(data1_id);
         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > d2 =
                                                  p->getPatchData(data2_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!d1.isNull());
#endif
            box = d1->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }
         
         dprod += d_patch_ops.dot(d1, d2, box, cv);
      }
   }

   dcomplex global_dot = tbox::MPI::sumReduction(dprod);
   return( global_dot );
}

template<int DIM> dcomplex HierarchyFaceDataOpsComplex<DIM>::integral(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   dcomplex local_integral = dcomplex(0.0,0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,dcomplex> > data =
                                                  p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > vol = p->getPatchData(vol_id);

#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!data.isNull());
         assert(!vol.isNull());
#endif

         hier::Box<DIM> box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   dcomplex global_integral = tbox::MPI::sumReduction(local_integral);
   return( global_integral );
}

}
}
#endif
