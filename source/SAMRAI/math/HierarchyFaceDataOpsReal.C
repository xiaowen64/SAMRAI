/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real face data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyFaceDataOpsReal_C
#define included_math_HierarchyFaceDataOpsReal_C

#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

template<class TYPE>
HierarchyFaceDataOpsReal<TYPE>::HierarchyFaceDataOpsReal(
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsReal<TYPE>()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
#endif
   d_hierarchy = hierarchy;
   if ((coarsest_level < 0) || (finest_level < 0)) {
      if (d_hierarchy->getNumberOfLevels() == 0) {
         d_coarsest_level = coarsest_level;
         d_finest_level = finest_level;
      } else {
         resetLevels(0, d_hierarchy->getFinestLevelNumber());
      }
   } else {
      resetLevels(coarsest_level, finest_level);
   }
}

template<class TYPE>
HierarchyFaceDataOpsReal<TYPE>::~HierarchyFaceDataOpsReal()
{
}

/*
 *************************************************************************
 *
 * Routines to set the hierarchy and level informtation.
 *
 *************************************************************************
 */

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::setPatchHierarchy(
   tbox::Pointer<hier::PatchHierarchy> hierarchy)
{
   TBOX_ASSERT(!hierarchy.isNull());

   d_hierarchy = hierarchy;
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   const tbox::Dimension& dim(d_hierarchy->getDim());

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   for (int d = 0; d < dim.getValue(); d++) {
      d_nonoverlapping_face_boxes[d].resizeArray(d_finest_level + 1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::BoxList face_boxes;

      for (int nd = 0; nd < dim.getValue(); nd++) {
         face_boxes = level->getBoxes();
         for (hier::BoxList::Iterator i(face_boxes); i != face_boxes.end(); ++i) {
            *i = pdat::FaceGeometry::toFaceBox(*i, nd);
         }
         hier::BoxUtilities::makeNonOverlappingBoxLists(
            d_nonoverlapping_face_boxes[nd][ln],
            face_boxes);
      }
   }
}

template<class TYPE>
const tbox::Pointer<hier::PatchHierarchy>
HierarchyFaceDataOpsReal<TYPE>::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *
 * The following are private and cannot be used, but they are defined
 * here for compilers that require that every template declaration have
 * a definition (a stupid requirement, if you ask me).
 *
 *************************************************************************
 */

template<class TYPE>
HierarchyFaceDataOpsReal<TYPE>::HierarchyFaceDataOpsReal(
   const HierarchyFaceDataOpsReal<TYPE>& foo):
   HierarchyDataOpsReal<TYPE>()
{
   NULL_USE(foo);
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::operator = (
   const HierarchyFaceDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Basic generic operations.
 *
 *************************************************************************
 */

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::copyData(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.copyData(dst, src, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer<pdat::FaceDataFactory<TYPE> >
   d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   TBOX_ASSERT(!d1fact.isNull());
   tbox::Pointer<pdat::FaceDataFactory<TYPE> >
   d2fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id);
   TBOX_ASSERT(!d2fact.isNull());
   TBOX_ASSERT(d1fact->getDepth() == d2fact->getDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   s << "Patch descriptor id = " << data_id << std::endl;
   s << "Factory = " << typeid(*d_hierarchy->getPatchDescriptor()->
                               getPatchDataFactory(data_id)).name()
     << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << std::endl;
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::setToScalar(
   const int data_id,
   const TYPE& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setToScalar(d, alpha, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Basic generic arithmetic operations.
 *
 *************************************************************************
 */

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::scale(
   const int dst_id,
   const TYPE& alpha,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::addScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.add(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.subtract(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.multiply(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.divide(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.reciprocal(dst, src, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::linearSum(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const TYPE& beta,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.linearSum(dst, alpha, src1, beta, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::axpy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.axpy(dst, alpha, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::axmy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.axmy(dst, alpha, src1, src2, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.abs(dst, src, box);
      }
   }
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::setRandomValues(
   const int data_id,
   const TYPE& width,
   const TYPE& low,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!data.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : data->getGhostBox());

         d_patch_ops.setRandomValues(data, width, low, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Generic norm and order operations.
 *
 *************************************************************************
 */

template<class TYPE>
int HierarchyFaceDataOpsReal<TYPE>::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   const tbox::Dimension& dim(d_hierarchy->getDim());

   int entries = 0;

   if (interior_only) {

      tbox::Pointer<pdat::FaceDataFactory<TYPE> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < dim.getValue(); nd++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_face_boxes[nd][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {
            for (int eb = 0; eb < dim.getValue(); eb++) {
               hier::BoxList::ConstIterator lb =
                  ((d_nonoverlapping_face_boxes[eb][ln])[il]).begin();
               for ( ; lb != ((d_nonoverlapping_face_boxes[eb][ln])[il]).end();
                    ++lb) {
                  entries += lb().size();
               }
            }
         }
      }

      entries *= dfact->getDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
         for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
            tbox::Pointer<pdat::FaceData<TYPE> > d =
               (*ip)->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      int global_entries = entries;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&entries, &global_entries, 1, MPI_INT, MPI_SUM);
      }
      entries = global_entries;

   }

   return entries;
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(vol_id >= 0);
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer<pdat::FaceData<double> > cv = p->getPatchData(vol_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!cv.isNull());
#endif
         hier::Box box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(data, cv, box);
      }
   }

   double global_sum = sum;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_sum;
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::L1Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         norm += d_patch_ops.L1Norm(data, box, cv);
      }
   }

   if (!local_only) {
      double global_norm = norm;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
      }
      norm = global_norm;
   }
   return norm;
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::L2Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   double norm_squared = HierarchyFaceDataOpsReal<TYPE>::dot(data_id,
         data_id,
         vol_id,
         local_only);

   return sqrt(norm_squared);
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::weightedL2Norm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm_squared = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer<pdat::FaceData<TYPE> > weight = p->getPatchData(wgt_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         double pnorm = d_patch_ops.weightedL2Norm(data, weight, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = norm_squared;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm_squared, &global_norm_squared, 1, MPI_DOUBLE, MPI_SUM);
   }
   return sqrt(global_norm_squared);
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::weightedRMSNorm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{

   double l2_norm = weightedL2Norm(data_id, wgt_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

template<class TYPE>
double HierarchyFaceDataOpsReal<TYPE>::maxNorm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data = p->getPatchData(data_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            box = data->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         norm = tbox::MathUtilities<double>::Max(norm,
               d_patch_ops.maxNorm(data, box, cv));
      }
   }

   if (!local_only) {
      double global_norm = norm;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX);
      }
      norm = global_norm;
   }
   return norm;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE dprod = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data1.isNull());
#endif
            box = data1->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         dprod += d_patch_ops.dot(data1, data2, box, cv);
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&dprod, 1, MPI_SUM);
      }
   }
   return dprod;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::integral(
   const int data_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE local_integral = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data =
            p->getPatchData(data_id);
         tbox::Pointer<pdat::FaceData<double> > vol = p->getPatchData(vol_id);

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!data.isNull());
         TBOX_ASSERT(!vol.isNull());
#endif

         hier::Box box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   TYPE global_integral = local_integral;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_integral, 1, MPI_SUM);
   }
   return global_integral;
}

/*
 *************************************************************************
 *
 * Generic miscellaneous operations for real data.
 *
 *************************************************************************
 */

template<class TYPE>
int HierarchyFaceDataOpsReal<TYPE>::computeConstrProdPos(
   const int data1_id,
   const int data2_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > data1 = p->getPatchData(data1_id);
         tbox::Pointer<pdat::FaceData<TYPE> > data2 = p->getPatchData(data2_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data1.isNull());
#endif
            box = data1->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         test = tbox::MathUtilities<int>::Min(test,
               d_patch_ops.computeConstrProdPos(data1, data2, box, cv));
      }
   }

   int global_test = test;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&test, &global_test, 1, MPI_INT, MPI_MIN);
   }
   return global_test;
}

template<class TYPE>
void HierarchyFaceDataOpsReal<TYPE>::compareToScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         d_patch_ops.compareToScalar(dst, src, alpha, box, cv);
      }
   }
}

template<class TYPE>
int HierarchyFaceDataOpsReal<TYPE>::testReciprocal(
   const int dst_id,
   const int src_id,
   const int vol_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::FaceData<TYPE> > src = p->getPatchData(src_id);
         tbox::Pointer<pdat::FaceData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            box = dst->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         test = tbox::MathUtilities<int>::Min(test,
               d_patch_ops.testReciprocal(dst, src, box, cv));
      }
   }

   int global_test = test;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&test, &global_test, 1, MPI_INT, MPI_MIN);
   }
   return global_test;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::maxPointwiseDivide(
   const int numer_id,
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE max = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > numer = p->getPatchData(numer_id);
         tbox::Pointer<pdat::FaceData<TYPE> > denom = p->getPatchData(denom_id);

         hier::Box box = p->getBox();

         max = tbox::MathUtilities<TYPE>::Max(max,
               d_patch_ops.maxPointwiseDivide(numer, denom, box));
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&max, 1, MPI_MAX);
      }
   }
   return max;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::minPointwiseDivide(
   const int numer_id,
   const int denom_id,
   bool local_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE min = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > numer = p->getPatchData(numer_id);
         tbox::Pointer<pdat::FaceData<TYPE> > denom = p->getPatchData(denom_id);

         hier::Box box = p->getBox();

         min = tbox::MathUtilities<TYPE>::Min(min,
               d_patch_ops.minPointwiseDivide(numer, denom, box));
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&min, 1, MPI_MIN);
      }
   }
   return min;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         minval = tbox::MathUtilities<TYPE>::Min(minval, d_patch_ops.min(d, box));
      }
   }

   TYPE global_min = minval;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_min, 1, MPI_MIN);
   }
   return global_min;
}

template<class TYPE>
TYPE HierarchyFaceDataOpsReal<TYPE>::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::FaceData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         maxval = tbox::MathUtilities<TYPE>::Max(maxval, d_patch_ops.max(d, box));
      }
   }

   TYPE global_max = maxval;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_max, 1, MPI_MAX);
   }
   return global_max;
}

}
}
#endif
