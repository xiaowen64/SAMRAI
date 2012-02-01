/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer side data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchySideDataOpsInteger_C
#define included_math_HierarchySideDataOpsInteger_C

#include "SAMRAI/math/HierarchySideDataOpsInteger.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

HierarchySideDataOpsInteger::HierarchySideDataOpsInteger(
   boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsInteger(),
   d_hierarchy(hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
#endif
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

HierarchySideDataOpsInteger::~HierarchySideDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Rotuines to set the hierarchy and level information.
 *
 *************************************************************************
 */

void HierarchySideDataOpsInteger::setPatchHierarchy(
   boost::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   TBOX_ASSERT(hierarchy);

   d_hierarchy = hierarchy;
}

void HierarchySideDataOpsInteger::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   const tbox::Dimension& dim(d_hierarchy->getDim());

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   for (int d = 0; d < dim.getValue(); d++) {
      d_nonoverlapping_side_boxes[d].resizeArray(d_finest_level + 1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      hier::BoxContainer side_boxes;

      for (int nd = 0; nd < dim.getValue(); nd++) {
         side_boxes = level->getBoxes();
         for (hier::BoxContainer::Iterator i(side_boxes); i != side_boxes.end(); ++i) {
            *i = pdat::SideGeometry::toSideBox(*i, nd);
         }
         hier::BoxUtilities::makeNonOverlappingBoxContainers(
            d_nonoverlapping_side_boxes[nd][ln],
            side_boxes);
      }
   }
}

const boost::shared_ptr<hier::PatchHierarchy>
HierarchySideDataOpsInteger::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *
 * Basic generic operations.
 *
 *************************************************************************
 */

int HierarchySideDataOpsInteger::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
   const tbox::Dimension& dim(d_hierarchy->getDim());

   int entries = 0;

   if (interior_only) {

      boost::shared_ptr<pdat::SideDataFactory<int> > dfact(
         d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id),
         boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(dfact);
#endif
      const hier::IntVector& directions = dfact->getDirectionVector();

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int dc = 0; dc < dim.getValue(); dc++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_side_boxes[dc][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {
            for (int eb = 0; eb < dim.getValue(); eb++) {
               if (directions(eb)) {
                  hier::BoxContainer::ConstIterator lb =
                     ((d_nonoverlapping_side_boxes[eb][ln])[il]).begin();
                  for ( ; lb != ((d_nonoverlapping_side_boxes[eb][ln])[il]).end();
                       ++lb) {
                     entries += lb().size();
                  }
               }
            }
         }
      }

      entries *= dfact->getDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
            boost::shared_ptr<pdat::SideData<int> > d(
               (*ip)->getPatchData(data_id),
               boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d);
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

void HierarchySideDataOpsInteger::copyData(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.copyData(d, s, box);
      }
   }
}

void HierarchySideDataOpsInteger::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   boost::shared_ptr<pdat::SideDataFactory<int> > d1fact(
      d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(d1fact);
   boost::shared_ptr<pdat::SideDataFactory<int> > d2fact(
      d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id),
   boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(d2fact);
   TBOX_ASSERT(d1fact->getDepth() == d2fact->getDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() ==
      d2fact->getGhostCellWidth());
   TBOX_ASSERT(d1fact->getDirectionVector() == d2fact->getDirectionVector());
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

void HierarchySideDataOpsInteger::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
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
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

void HierarchySideDataOpsInteger::setToScalar(
   const int data_id,
   const int& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
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

void HierarchySideDataOpsInteger::scale(
   const int dst_id,
   const int& alpha,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > dst(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(dst);
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

void HierarchySideDataOpsInteger::addScalar(
   const int dst_id,
   const int src_id,
   const int& alpha,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > dst(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(dst);
#endif
         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

void HierarchySideDataOpsInteger::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

void HierarchySideDataOpsInteger::linearSum(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int& beta,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::axpy(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::axmy(
   const int dst_id,
   const int& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

void HierarchySideDataOpsInteger::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<int> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.abs(d, src, box);
      }
   }
}

int HierarchySideDataOpsInteger::min(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int minval = tbox::MathUtilities<int>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
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

int HierarchySideDataOpsInteger::max(
   const int data_id,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int maxval = -(tbox::MathUtilities<int>::getMax());

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
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

void HierarchySideDataOpsInteger::setRandomValues(
   const int data_id,
   const int& width,
   const int& low,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::SideData<int> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d);
#endif
         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setRandomValues(d, width, low, box);
      }
   }
}

}
}
#endif
