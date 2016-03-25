//
// File:	HierarchyFaceDataOpsReal.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Templated operations for real face data on multiple levels.
//

#ifndef included_math_HierarchyFaceDataOpsReal_C
#define included_math_HierarchyFaceDataOpsReal_C

#include "HierarchyFaceDataOpsReal.h"
#include "BoxUtilities.h"
#include "PatchDescriptor.h"
#include "FaceDataFactory.h"
#include "tbox/IEEE.h"
#include "tbox/MPI.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
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

template<int DIM, class TYPE>
HierarchyFaceDataOpsReal<DIM,TYPE>::HierarchyFaceDataOpsReal(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
: HierarchyDataOpsReal<DIM,TYPE>()
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
 
template<int DIM, class TYPE>
HierarchyFaceDataOpsReal<DIM,TYPE>::~HierarchyFaceDataOpsReal()
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
void HierarchyFaceDataOpsReal<DIM,TYPE>::setPatchHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::resetLevels(
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

template<int DIM, class TYPE>
const tbox::Pointer< hier::PatchHierarchy<DIM> >
HierarchyFaceDataOpsReal<DIM,TYPE>::getPatchHierarchy() const
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
HierarchyFaceDataOpsReal<DIM,TYPE>::HierarchyFaceDataOpsReal(
   const HierarchyFaceDataOpsReal<DIM,TYPE>& foo)
:  HierarchyDataOpsReal<DIM,TYPE>()
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::operator=(
   const HierarchyFaceDataOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

/*
*************************************************************************
*                                                                       *
* Basic generic operations.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::copyData(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull()); 
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.copyData(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer< pdat::FaceDataFactory<DIM,TYPE> >
      d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   assert(!d1fact.isNull());
   tbox::Pointer< pdat::FaceDataFactory<DIM,TYPE> >
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

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::printData(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::setToScalar(
   const int data_id,
   const TYPE& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > d = p->getPatchData(data_id);
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

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::scale(
   const int dst_id,
   const TYPE& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::addScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::add(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.add(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::subtract(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.subtract(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::multiply(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.multiply(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::divide(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.divide(dst, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::reciprocal(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.reciprocal(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::linearSum(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const TYPE& beta,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.linearSum(dst, alpha, src1, beta, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::axpy(
   const int dst_id,
   const TYPE& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );
         
         d_patch_ops.axpy(dst, alpha, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::axmy(
   const int dst_id,
   const TYPE& alpha,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.axmy(dst, alpha, src1, src2, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::abs(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!dst.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : dst->getGhostBox() );

         d_patch_ops.abs(dst, src, box);
      }
   }
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::setRandomValues(
   const int data_id,
   const TYPE& width,
   const TYPE& low,
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!data.isNull());
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
int HierarchyFaceDataOpsReal<DIM,TYPE>::numberOfEntries(
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

      tbox::Pointer< pdat::FaceDataFactory<DIM,TYPE> >
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
            tbox::Pointer< pdat::FaceData<DIM,TYPE> > d = 
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

template<int DIM, class TYPE>
double HierarchyFaceDataOpsReal<DIM,TYPE>::sumControlVolumes(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data  = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv  = p->getPatchData(vol_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!cv.isNull());
#endif
         hier::Box<DIM> box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(data, cv, box);
      }
   }

   double global_sum = tbox::MPI::sumReduction(sum);
   return( global_sum );
}

template<int DIM, class TYPE>
double HierarchyFaceDataOpsReal<DIM,TYPE>::L1Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm += d_patch_ops.L1Norm(data, box, cv);
      }
   }

   if ( !local_only ) {
      norm = tbox::MPI::sumReduction(norm);
   }
   return( norm );
}

template<int DIM, class TYPE>
double HierarchyFaceDataOpsReal<DIM,TYPE>::L2Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   double norm_squared = HierarchyFaceDataOpsReal<DIM,TYPE>::dot(data_id,
                                                                   data_id,
                                                                   vol_id,
                                                                   local_only);

   return( sqrt(norm_squared) );
}

template<int DIM, class TYPE>
double HierarchyFaceDataOpsReal<DIM,TYPE>::weightedL2Norm(
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data   = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > weight = p->getPatchData(wgt_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         double pnorm = d_patch_ops.weightedL2Norm(data, weight, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = tbox::MPI::sumReduction(norm_squared);
   return( sqrt(global_norm_squared) );
}

template<int DIM, class TYPE>
double HierarchyFaceDataOpsReal<DIM,TYPE>::RMSNorm(
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
double HierarchyFaceDataOpsReal<DIM,TYPE>::weightedRMSNorm(
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
double HierarchyFaceDataOpsReal<DIM,TYPE>::maxNorm(
   const int data_id,
   const int vol_id,
   bool local_only) const
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!data.isNull());
#endif
            box = data->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         norm = tbox::Utilities::dmax(norm, d_patch_ops.maxNorm(data, box, cv));
      }
   }

   if ( !local_only ) {
      norm = tbox::MPI::maxReduction(norm);
   }
   return( norm );
}

template<int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE dprod = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!data1.isNull());
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
      TYPE global_dot = tbox::MPI::sumReduction(dprod);
      return( global_dot );
   }
}

template<int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::integral(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE local_integral = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data =
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

   TYPE global_integral = tbox::MPI::sumReduction(local_integral);
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
int HierarchyFaceDataOpsReal<DIM,TYPE>::computeConstrProdPos(
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

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!data1.isNull());
#endif
            box = data1->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         test = tbox::Utilities::imin( test, 
                   d_patch_ops.computeConstrProdPos(data1, data2, box, cv) );
      }
   }

   int global_test = tbox::MPI::minReduction(test);
   return( global_test );
}

template<int DIM, class TYPE>
void HierarchyFaceDataOpsReal<DIM,TYPE>::compareToScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const int vol_id) const
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

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src  = p->getPatchData(src_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }
  
         d_patch_ops.compareToScalar(dst, src, alpha, box, cv);
      }
   }
}

template<int DIM, class TYPE>
int HierarchyFaceDataOpsReal<DIM,TYPE>::testReciprocal(
   const int dst_id, 
   const int src_id, 
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > dst  = p->getPatchData(dst_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > src  = p->getPatchData(src_id);
         tbox::Pointer< pdat::FaceData<DIM,double> > cv;

         hier::Box<DIM> box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv  = p->getPatchData(vol_id);
         }

         test = tbox::Utilities::imin(test,
                   d_patch_ops.testReciprocal(dst, src, box, cv) );
      }
   }

   int global_test = tbox::MPI::minReduction(test);
   return( global_test );
}

template<int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const int numer_id, 
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE max = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > numer  = p->getPatchData(numer_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > denom  = p->getPatchData(denom_id);

         hier::Box<DIM> box = p->getBox();

         max = tbox::Utilities::dmax(max,
                   d_patch_ops.maxPointwiseDivide(numer, denom, box));
      }
   }

   if ( !local_only ) {
      max = tbox::MPI::maxReduction(max);
   }
   return( max );
}

template <int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::minPointwiseDivide(
   const int numer_id, 
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE min = tbox::IEEE::getFLT_MAX();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > numer  = p->getPatchData(numer_id);
         tbox::Pointer< pdat::FaceData<DIM,TYPE> > denom  = p->getPatchData(denom_id);

         hier::Box<DIM> box = p->getBox();

         min = tbox::Utilities::dmin(min,
                   d_patch_ops.minPointwiseDivide(numer, denom, box));
      }
   }

   if ( !local_only ) {
      min = tbox::MPI::minReduction(min);
   }
   return( min );
}

template<int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         minval = tbox::MathUtilities<TYPE>::Min( minval, d_patch_ops.min(d, box) );
      }
   }

   TYPE global_min = tbox::MPI::minReduction(minval);
   return( global_min );
}

template<int DIM, class TYPE>
TYPE HierarchyFaceDataOpsReal<DIM,TYPE>::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > p = level->getPatch(ip());

         tbox::Pointer< pdat::FaceData<DIM,TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d.isNull());
#endif
         hier::Box<DIM> box = ( interior_only ? p->getBox() : d->getGhostBox() );

         maxval = tbox::MathUtilities<TYPE>::Max( maxval, d_patch_ops.max(d, box) );
      }
   }

   TYPE global_max = tbox::MPI::maxReduction(maxval);
   return( global_max );
}

}
}
#endif
