/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex node data on multiple levels. 
 *
 ************************************************************************/

#ifndef included_math_HierarchyNodeDataOpsComplex_C
#define included_math_HierarchyNodeDataOpsComplex_C

#include "SAMRAI/math/HierarchyNodeDataOpsComplex.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

HierarchyNodeDataOpsComplex::HierarchyNodeDataOpsComplex(
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsComplex()
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

HierarchyNodeDataOpsComplex::~HierarchyNodeDataOpsComplex()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Routines to set the hierarchy and level information.                  *
 *                                                                       *
 *************************************************************************
 */

void HierarchyNodeDataOpsComplex::setPatchHierarchy(
   tbox::Pointer<hier::PatchHierarchy> hierarchy)
{
   TBOX_ASSERT(!hierarchy.isNull());

   d_hierarchy = hierarchy;
}

void HierarchyNodeDataOpsComplex::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   d_nonoverlapping_node_boxes.resizeArray(d_finest_level + 1);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::BoxList node_boxes = level->getBoxes();

      for (hier::BoxList::Iterator i(node_boxes); i; i++) {
         *i = pdat::NodeGeometry::toNodeBox(*i);
      }
      hier::BoxUtilities::makeNonOverlappingBoxLists(
         d_nonoverlapping_node_boxes[ln],
         node_boxes);
   }
}

const tbox::Pointer<hier::PatchHierarchy>
HierarchyNodeDataOpsComplex::getPatchHierarchy() const
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

void HierarchyNodeDataOpsComplex::copyData(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.copyData(d, s, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Pointer<pdat::NodeDataFactory<dcomplex> >
   d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
   TBOX_ASSERT(!d1fact.isNull());
   tbox::Pointer<pdat::NodeDataFactory<dcomplex> >
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

void HierarchyNodeDataOpsComplex::printData(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

void HierarchyNodeDataOpsComplex::setToScalar(
   const int data_id,
   const dcomplex& alpha,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
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

void HierarchyNodeDataOpsComplex::scale(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::addScalar(
   const int dst_id,
   const int src_id,
   const dcomplex& alpha,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > dst = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dst.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::add(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::subtract(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::multiply(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::divide(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::reciprocal(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::linearSum(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const dcomplex& beta,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::axpy(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::axmy(
   const int dst_id,
   const dcomplex& alpha,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s1 = p->getPatchData(src1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > s2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::abs(
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

         tbox::Pointer<pdat::NodeData<double> > d = p->getPatchData(dst_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.abs(d, src, box);
      }
   }
}

void HierarchyNodeDataOpsComplex::setRandomValues(
   const int data_id,
   const dcomplex& width,
   const dcomplex& low,
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!d.isNull());
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

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

int HierarchyNodeDataOpsComplex::numberOfEntries(
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

   int entries = 0;

   if (interior_only) {

      tbox::Pointer<pdat::NodeDataFactory<dcomplex> >
      dfact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!dfact.isNull());
#endif

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONa
         TBOX_ASSERT(npatches == d_nonoverlapping_node_boxes[ln].getSize());
#endif
         for (int il = 0; il < npatches; il++) {
            tbox::List<hier::Box>::Iterator lb =
               ((d_nonoverlapping_node_boxes[ln])[il]).listStart();
            for ( ; lb; lb++) {
               entries += lb().size();
            }
         }
      }

      entries *= dfact->getDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
         for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
            tbox::Pointer<pdat::NodeData<dcomplex> > d =
               (*ip)->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      int global_entries = entries;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&entries, &global_entries, 1, MPI_INT, MPI_MAX);
      }
      entries = global_entries;

   }

   return entries;
}

double HierarchyNodeDataOpsComplex::sumControlVolumes(
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

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer<pdat::NodeData<double> > cv = p->getPatchData(vol_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!cv.isNull());
#endif
         hier::Box box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(d, cv, box);
      }
   }

   double global_sum = sum;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_sum;
}

double HierarchyNodeDataOpsComplex::L1Norm(
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

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer<pdat::NodeData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            box = d->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         norm += d_patch_ops.L1Norm(d, box, cv);
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_norm;
}

double HierarchyNodeDataOpsComplex::L2Norm(
   const int data_id,
   const int vol_id) const
{
   dcomplex dotprod = HierarchyNodeDataOpsComplex::dot(data_id,
         data_id,
         vol_id);

   return sqrt(real(dotprod));
}

double HierarchyNodeDataOpsComplex::weightedL2Norm(
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

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > w = p->getPatchData(wgt_id);
         tbox::Pointer<pdat::NodeData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            box = d->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         double pnorm = d_patch_ops.weightedL2Norm(d, w, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = norm_squared;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm_squared, &global_norm_squared, 1, MPI_DOUBLE, MPI_SUM);
   }
   return sqrt(global_norm_squared);
}

double HierarchyNodeDataOpsComplex::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

double HierarchyNodeDataOpsComplex::weightedRMSNorm(
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

double HierarchyNodeDataOpsComplex::maxNorm(
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

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::NodeData<dcomplex> > d = p->getPatchData(data_id);
         tbox::Pointer<pdat::NodeData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            box = d->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         norm = tbox::MathUtilities<double>::Max(norm,
               d_patch_ops.maxNorm(d, box, cv));
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_norm;
}

dcomplex HierarchyNodeDataOpsComplex::dot(
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

   dcomplex dprod = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::NodeData<dcomplex> > d1 =
            p->getPatchData(data1_id);
         tbox::Pointer<pdat::NodeData<dcomplex> > d2 =
            p->getPatchData(data2_id);
         tbox::Pointer<pdat::NodeData<double> > cv;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d1.isNull());
#endif
            box = d1->getGhostBox();
            cv = p->getPatchData(vol_id);
         }

         dprod += d_patch_ops.dot(d1, d2, box, cv);
      }
   }

   dcomplex global_dot = dprod;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&dprod, &global_dot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM);
   }
   return global_dot;
}

dcomplex HierarchyNodeDataOpsComplex::integral(
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

   dcomplex local_integral = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch> p = *ip;

         tbox::Pointer<pdat::NodeData<dcomplex> > data =
            p->getPatchData(data_id);
         tbox::Pointer<pdat::NodeData<double> > vol = p->getPatchData(vol_id);

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!data.isNull());
         TBOX_ASSERT(!vol.isNull());
#endif

         hier::Box box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   dcomplex global_integral = local_integral;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&local_integral, &global_integral, 1, MPI_DOUBLE_COMPLEX, MPI_SUM);
   }
   return global_integral;
}

}
}
#endif
