/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer edge data on multiple levels. 
 *
 ************************************************************************/

#ifndef included_math_HierarchyEdgeDataOpsInteger_C
#define included_math_HierarchyEdgeDataOpsInteger_C

#include "SAMRAI/math/HierarchyEdgeDataOpsInteger.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

HierarchyEdgeDataOpsInteger::HierarchyEdgeDataOpsInteger(
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsInteger()
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

HierarchyEdgeDataOpsInteger::~HierarchyEdgeDataOpsInteger()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Rotuines to set the hierarchy and level information.                  *
 *                                                                       *
 *************************************************************************
 */

void HierarchyEdgeDataOpsInteger::setPatchHierarchy(
   tbox::Pointer<hier::PatchHierarchy> hierarchy)
{
   TBOX_ASSERT(!hierarchy.isNull());

   d_hierarchy = hierarchy;
}

void HierarchyEdgeDataOpsInteger::resetLevels(
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
      d_nonoverlapping_edge_boxes[d].resizeArray(d_finest_level + 1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::BoxList edge_boxes(dim);

      for (int nd = 0; nd < dim.getValue(); nd++) {
         edge_boxes = level->getBoxes();
         for (hier::BoxList::Iterator i(edge_boxes); i; i++) {
            *i = pdat::EdgeGeometry::toEdgeBox(*i, nd);
         }
         hier::BoxUtilities::makeNonOverlappingBoxLists(
            d_nonoverlapping_edge_boxes[nd][ln],
            edge_boxes);
      }
   }
}

const tbox::Pointer<hier::PatchHierarchy>
HierarchyEdgeDataOpsInteger::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *                                                                       *
 * Basic generic operations.                                             *
 *                                                                       *
 *************************************************************************
 */

int HierarchyEdgeDataOpsInteger::numberOfEntries(
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

      tbox::Pointer<pdat::EdgeDataFactory<int> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < dim.getValue(); nd++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_edge_boxes[nd][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {
            tbox::List<hier::Box>::Iterator lb;

            for (int eb = 0; eb < dim.getValue(); eb++) {
               lb = ((d_nonoverlapping_edge_boxes[eb][ln])[il]).listStart();
               for ( ; lb; lb++) {
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
            tbox::Pointer<pdat::EdgeData<int> > d =
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

void HierarchyEdgeDataOpsInteger::copyData(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.copyData(d, s, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer<pdat::EdgeDataFactory<int> >
   d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   TBOX_ASSERT(!d1fact.isNull());
   tbox::Pointer<pdat::EdgeDataFactory<int> >
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

void HierarchyEdgeDataOpsInteger::printData(
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
   s << "Factory = " << typeid(d_hierarchy->getPatchDescriptor()->
                               getPatchDataFactory(data_id)).name()
     << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << std::endl;
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

void HierarchyEdgeDataOpsInteger::setToScalar(
   const int data_id,
   const int& alpha,
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(data_id);
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
 *                                                                       *
 * Basic generic arithmetic operations.                                  *
 *                                                                       *
 *************************************************************************
 */

void HierarchyEdgeDataOpsInteger::scale(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer<pdat::EdgeData<int> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::addScalar(
   const int dst_id,
   const int src_id,
   const int& alpha,
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

         tbox::Pointer<pdat::EdgeData<int> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::add(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::subtract(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::multiply(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::divide(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::reciprocal(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::linearSum(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int& beta,
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::axpy(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::axmy(
   const int dst_id,
   const int& alpha,
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::EdgeData<int> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchyEdgeDataOpsInteger::abs(
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::EdgeData<int> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.abs(d, src, box);
      }
   }
}

int HierarchyEdgeDataOpsInteger::min(
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

   int minval = tbox::MathUtilities<int>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         minval = tbox::MathUtilities<int>::Min(minval,
               d_patch_ops.min(d, box));
      }
   }

   int global_min = minval;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&minval, &global_min, 1, MPI_INT, MPI_MIN);
   }
   return global_min;
}

int HierarchyEdgeDataOpsInteger::max(
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

   int maxval = -(tbox::MathUtilities<int>::getMax());

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         maxval = tbox::MathUtilities<int>::Max(maxval,
               d_patch_ops.min(d, box));
      }
   }

   int global_max = maxval;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&maxval, &global_max, 1, MPI_INT, MPI_MAX);
   }
   return global_max;
}

void HierarchyEdgeDataOpsInteger::setRandomValues(
   const int data_id,
   const int& width,
   const int& low,
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

         tbox::Pointer<pdat::EdgeData<int> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setRandomValues(d, width, low, box);
      }
   }
}

}
}
#endif
