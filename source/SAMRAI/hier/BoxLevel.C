/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of boxes in a box_level of a distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevel_C
#define included_hier_BoxLevel_C

#include "SAMRAI/hier/BoxLevel.h"

#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxLevel.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int BoxLevel::HIER_MAPPED_BOX_LEVEL_VERSION = 0;
const int BoxLevel::MAPPED_BOX_LEVEL_NUMBER_OF_STATS = 20;

tbox::Pointer<tbox::Timer> BoxLevel::t_initialize_private;
tbox::Pointer<tbox::Timer> BoxLevel::t_acquire_remote_boxes;
tbox::Pointer<tbox::Timer> BoxLevel::t_cache_global_reduced_data;

const LocalId BoxLevel::s_negative_one_local_id(-1);

tbox::StartupShutdownManager::Handler
BoxLevel::s_initialize_finalize_handler(
   BoxLevel::initializeCallback,
   0,
   0,
   BoxLevel::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

BoxLevel::BoxLevel():
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_ratio(tbox::Dimension::getInvalidDimension(), 0),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(-1),

   d_max_number_of_boxes(-1),
   d_min_number_of_boxes(-1),
   d_max_number_of_cells(-1),
   d_min_number_of_cells(-1),

   d_local_max_box_size(0),
   d_global_max_box_size(0),
   d_local_min_box_size(0),
   d_global_min_box_size(0),

   d_local_bounding_box(0),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(0),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(NULL),
   d_persistent_overlap_connectors(NULL),
   d_handle(NULL),
   d_grid_geometry(tbox::ConstPointer<GridGeometry>(NULL))
{
   // This ctor should never be invoked.
   TBOX_ERROR("Somehow, we entered code that was never meant to be used.");
}

BoxLevel::BoxLevel(
   const tbox::Dimension& dim):

   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_ratio(dim, 0),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(-1),

   d_max_number_of_boxes(-1),
   d_min_number_of_boxes(-1),
   d_max_number_of_cells(-1),
   d_min_number_of_cells(-1),

   d_local_max_box_size(0),
   d_global_max_box_size(0),
   d_local_min_box_size(0),
   d_global_min_box_size(0),

   d_local_bounding_box(0),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(0),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(NULL),
   d_persistent_overlap_connectors(NULL),
   d_handle(NULL),
   d_grid_geometry(tbox::ConstPointer<GridGeometry>(NULL))
{
}

BoxLevel::BoxLevel(
   const BoxLevel& rhs):
   tbox::DescribedClass(),
   d_mpi(rhs.d_mpi),
   d_boxes(rhs.d_boxes),
   d_global_boxes(rhs.d_global_boxes),
   d_ratio(rhs.d_ratio),

   d_local_number_of_cells(rhs.d_local_number_of_cells),
   d_global_number_of_cells(rhs.d_global_number_of_cells),
   d_local_number_of_boxes(rhs.d_local_number_of_boxes),
   d_global_number_of_boxes(rhs.d_global_number_of_boxes),

   d_max_number_of_boxes(rhs.d_max_number_of_boxes),
   d_min_number_of_boxes(rhs.d_min_number_of_boxes),
   d_max_number_of_cells(rhs.d_max_number_of_cells),
   d_min_number_of_cells(rhs.d_min_number_of_cells),

   d_local_max_box_size(rhs.d_local_max_box_size),
   d_global_max_box_size(rhs.d_global_max_box_size),
   d_local_min_box_size(rhs.d_local_min_box_size),
   d_global_min_box_size(rhs.d_global_min_box_size),

   d_local_bounding_box(rhs.d_local_bounding_box),
   d_local_bounding_box_up_to_date(rhs.d_local_bounding_box_up_to_date),
   d_global_bounding_box(rhs.d_global_bounding_box),
   d_global_data_up_to_date(rhs.d_global_data_up_to_date),

   d_parallel_state(rhs.d_parallel_state),
   d_globalized_version(NULL),
   d_persistent_overlap_connectors(NULL),
   d_handle(NULL),
   d_grid_geometry(rhs.d_grid_geometry)
{
   // This cannot be the first constructor call, so no need to set timers.
}

BoxLevel::BoxLevel(
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state):
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_ratio(ratio),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(-1),

   d_max_number_of_boxes(-1),
   d_min_number_of_boxes(-1),
   d_max_number_of_cells(-1),
   d_min_number_of_cells(-1),

   d_local_max_box_size(0),
   d_global_max_box_size(0),
   d_local_min_box_size(0),
   d_global_min_box_size(0),

   d_local_bounding_box(0),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(0),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(NULL),
   d_persistent_overlap_connectors(NULL),
   d_handle(NULL),
   d_grid_geometry(tbox::ConstPointer<GridGeometry>(NULL))
{
   initialize(ratio, grid_geom, mpi, parallel_state);
}

BoxLevel::~BoxLevel()
{
   clear();
   if (d_persistent_overlap_connectors != NULL) {
      delete d_persistent_overlap_connectors;
      d_persistent_overlap_connectors = NULL;
   }
}

void BoxLevel::initialize(
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   d_boxes.clear();
   d_boxes.order();
   initializePrivate(
      ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void BoxLevel::swapInitialize(
   BoxContainer& boxes,
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   TBOX_ASSERT(&boxes != &d_boxes);   // Library error if this fails.
   d_boxes.swap(boxes);
   initializePrivate(ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void BoxLevel::finalize()
{

   // Erase non-local Boxes, if any, from d_boxes.
   for (BoxContainer::Iterator mbi = d_boxes.begin();
        mbi != d_boxes.end(); /* incremented in loop */) {
      if (mbi->getOwnerRank() != d_mpi.getRank()) {
         d_boxes.erase(mbi++);
      } else {
         ++mbi;
      }
   }

   computeLocalRedundantData();
   return;
}

void BoxLevel::initializePrivate(
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);
   t_initialize_private->start();

   clearForBoxChanges();

   d_mpi = mpi;
   d_grid_geometry = grid_geom;

   if (parallel_state == DISTRIBUTED) {
      d_global_boxes.clear();
   } else {
      d_global_boxes = d_boxes;
   }

   // Erase non-local Boxes, if any, from d_boxes.
   for (BoxContainer::Iterator mbi(d_boxes.begin());
        mbi != d_boxes.end(); /* incremented in loop */) {
      if (mbi->getOwnerRank() != d_mpi.getRank()) {
         d_boxes.erase(mbi++);
      } else {
         ++mbi;
      }
   }

   d_ratio = ratio;
   d_parallel_state = parallel_state;
   d_global_number_of_cells = -1;
   d_global_number_of_boxes = -1;
   d_local_bounding_box_up_to_date = false;
   d_global_data_up_to_date = false;
   computeLocalRedundantData();

   t_initialize_private->stop();
}

/*
 ***********************************************************************
 * Clear data and reset them to unusuable values.
 *
 * Note: don't use IntVector::getOne here, because SAMRAI may have
 * already shut down.
 ***********************************************************************
 */
void BoxLevel::removePeriodicImageBoxes()
{
   if (isInitialized()) {
      clearForBoxChanges();
      d_boxes.removePeriodicImageBoxes();
      if ( d_parallel_state == GLOBALIZED ) {
         d_global_boxes.removePeriodicImageBoxes();
      }
   }
}

/*
 ***********************************************************************
 * Clear data and reset them to unusuable values.
 *
 * Note: don't use IntVector::getOne here, because SAMRAI may have
 * already shut down.
 ***********************************************************************
 */
void BoxLevel::clear()
{
   if (isInitialized()) {
      clearForBoxChanges();
      d_mpi = tbox::SAMRAI_MPI(tbox::SAMRAI_MPI::commNull);
      d_boxes.clear();
      d_global_boxes.clear();
      d_ratio(0) = 0;
      d_local_number_of_cells = 0;
      d_global_number_of_cells = -1;
      d_local_number_of_boxes = 0;
      d_global_number_of_boxes = -1;
      d_local_bounding_box.clear();
      d_local_bounding_box_up_to_date = false;
      d_global_bounding_box.clear();
      d_global_data_up_to_date = false;
      d_local_max_box_size.clear();
      d_local_min_box_size.clear();
      d_global_max_box_size.clear();
      d_global_min_box_size.clear();
      d_parallel_state = DISTRIBUTED;
      d_grid_geometry = tbox::ConstPointer<GridGeometry>(NULL);
   }
}

void BoxLevel::swap(
   BoxLevel& level_a,
   BoxLevel& level_b)
{

   if (&level_a != &level_b) {
      if (level_a.isInitialized() && level_b.isInitialized()) {
         TBOX_DIM_ASSERT_CHECK_ARGS2(level_a, level_b);
      }

      level_a.clearPersistentOverlapConnectors();
      level_b.clearPersistentOverlapConnectors();

      level_a.detachMyHandle();
      level_b.detachMyHandle();

      // Swap objects supporting swap operation.
      level_a.d_boxes.swap(level_b.d_boxes);
      level_a.d_global_boxes.swap(level_b.d_global_boxes);
      level_a.d_local_bounding_box.swap(level_b.d_local_bounding_box);
      level_a.d_local_min_box_size.swap(level_b.d_local_min_box_size);
      level_a.d_local_max_box_size.swap(level_b.d_local_max_box_size);
      level_a.d_global_bounding_box.swap(level_b.d_global_bounding_box);

      // Swap objects not supporting swap operation.

      int tmpint;
      bool tmpbool;
      IntVector tmpvec(level_a.getDim());
      Box tmpbox(level_a.getDim());
      ParallelState tmpstate;
      const BoxLevel* tmpmbl;
      tbox::SAMRAI_MPI tmpmpi(tbox::SAMRAI_MPI::commNull);
      tbox::ConstPointer<GridGeometry> tmpgridgeom(level_a.getGridGeometry());

      tmpstate = level_a.d_parallel_state;
      level_a.d_parallel_state = level_b.d_parallel_state;
      level_b.d_parallel_state = tmpstate;

      tmpmpi = level_a.d_mpi;
      level_a.d_mpi = level_b.d_mpi;
      level_b.d_mpi = tmpmpi;

      tmpvec = level_a.d_ratio;
      level_a.d_ratio = level_b.d_ratio;
      level_b.d_ratio = tmpvec;

      tmpint = static_cast<int>(level_a.d_local_number_of_cells);
      level_a.d_local_number_of_cells = level_b.d_local_number_of_cells;
      level_b.d_local_number_of_cells = tmpint;

      tmpint = level_a.d_global_number_of_cells;
      level_a.d_global_number_of_cells = level_b.d_global_number_of_cells;
      level_b.d_global_number_of_cells = tmpint;

      tmpint = static_cast<int>(level_a.d_local_number_of_boxes);
      level_a.d_local_number_of_boxes = level_b.d_local_number_of_boxes;
      level_b.d_local_number_of_boxes = tmpint;

      tmpint = level_a.d_global_number_of_boxes;
      level_a.d_global_number_of_boxes = level_b.d_global_number_of_boxes;
      level_b.d_global_number_of_boxes = tmpint;

      tmpbool = level_a.d_local_bounding_box_up_to_date;
      level_a.d_local_bounding_box_up_to_date = level_b.d_local_bounding_box_up_to_date;
      level_b.d_local_bounding_box_up_to_date = tmpbool;

      tmpbool = level_a.d_global_data_up_to_date;
      level_a.d_global_data_up_to_date = level_b.d_global_data_up_to_date;
      level_b.d_global_data_up_to_date = tmpbool;

      tmpmbl = level_a.d_globalized_version;
      level_a.d_globalized_version = level_b.d_globalized_version;
      level_b.d_globalized_version = tmpmbl;

      tmpgridgeom = level_b.d_grid_geometry;
      level_a.d_grid_geometry = level_b.d_grid_geometry;
      level_b.d_grid_geometry = tmpgridgeom;
   }
}

void BoxLevel::computeLocalRedundantData()
{
   const IntVector max_vec(d_ratio.getDim(), tbox::MathUtilities<int>::getMax());
   const IntVector& zero_vec = IntVector::getZero(d_ratio.getDim());
   const int nblocks = d_grid_geometry->getNumberBlocks();

   d_local_number_of_boxes = 0;
   d_local_number_of_cells = 0;

   if (int(d_local_bounding_box.size()) != nblocks) {
      d_local_bounding_box.resize(nblocks, hier::Box(d_grid_geometry->getDim()));
      d_local_min_box_size.resize(nblocks, max_vec);
      d_local_max_box_size.resize(nblocks, zero_vec);
   }

   for (RealBoxConstIterator ni(d_boxes); ni.isValid(); ++ni) {

      int block_num = ni->getBlockId().getBlockValue();
      const IntVector boxdim(ni->numberCells());
      ++d_local_number_of_boxes;
      d_local_number_of_cells += boxdim.getProduct();
      d_local_bounding_box[block_num] += *ni;
      d_local_min_box_size[block_num].min(boxdim);
      d_local_max_box_size[block_num].max(boxdim);

   }

   d_local_bounding_box_up_to_date = true;
   d_global_data_up_to_date = false;
}

/*
 ****************************************************************************
 * Perform global reductions to get characteristics of the global data
 * without globalizing.  Data that can be reduced are combined into
 * arrays for reduction.  We do one sum reduction and one max reduction.
 * All the data we need fall into one or the other so it's all we need to
 * do.
 ****************************************************************************
 */
void BoxLevel::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isInitialized());

   if (d_global_data_up_to_date) {
      return;
   }

   t_cache_global_reduced_data->barrierAndStart();

   const int nblocks = d_grid_geometry->getNumberBlocks();

   /*
    * Sum reduction is used to compute the global sums of box count
    * and cell count.
    */
   if (d_parallel_state == GLOBALIZED) {
      d_global_number_of_boxes = 0;
      d_global_number_of_cells = 0;
      for (RealBoxConstIterator ni(d_global_boxes);
           ni.isValid();
           ++ni) {
         ++d_global_number_of_boxes;
         d_global_number_of_cells += ni->size();
      }
   } else {
      if (d_mpi.getSize() > 1) {
         int tmpa[2], tmpb[2];
         tmpa[0] = getLocalNumberOfBoxes();
         tmpa[1] = getLocalNumberOfCells();

         TBOX_ASSERT(tmpa[0] >= 0);
         TBOX_ASSERT(tmpa[1] >= 0);

         d_mpi.Allreduce(tmpa,
            tmpb,                        // Better to use MPI_IN_PLACE, but not some MPI's do not support.
            2,
            MPI_INT,
            MPI_SUM);
         d_global_number_of_boxes = tmpb[0];
         d_global_number_of_cells = tmpb[1];
      } else {
         d_global_number_of_boxes = getLocalNumberOfBoxes();
         d_global_number_of_cells = getLocalNumberOfCells();
      }

      TBOX_ASSERT(d_global_number_of_boxes >= 0);
      TBOX_ASSERT(d_global_number_of_cells >= 0);
   }

   if (int(d_global_bounding_box.size()) != nblocks) {
      d_global_bounding_box.resize(nblocks, hier::Box(getDim()));
      d_global_min_box_size.resize(nblocks, hier::IntVector(getDim()));
      d_global_max_box_size.resize(nblocks, hier::IntVector(getDim()));
   }

   /*
    * Max reduction is used to compute max/min box counts, max/min
    * cell counts, max/min box sizes, and bounding boxes.
    */
   if (d_mpi.getSize() == 1) {

      d_global_bounding_box = d_local_bounding_box;
      d_max_number_of_boxes = d_min_number_of_boxes = getLocalNumberOfBoxes();
      d_max_number_of_cells = d_min_number_of_cells = getLocalNumberOfCells();
      d_global_max_box_size = d_local_max_box_size;
      d_global_min_box_size = d_local_min_box_size;

   } else {

      if (d_mpi.getSize() > 1) {
         const tbox::Dimension& dim(getDim());

         std::vector<int> send_mesg;
         send_mesg.reserve(nblocks * 4 * dim.getValue() + 4);
         int bishift;
         for (int bn = 0; bn < nblocks; ++bn) {
            bishift = bn * nblocks;
            for (int i = 0; i < dim.getValue(); ++i) {
               send_mesg.push_back(-d_local_bounding_box[bn].lower()[i]);
               send_mesg.push_back(d_local_bounding_box[bn].upper()[i]);
               send_mesg.push_back(-d_local_min_box_size[bn][i]);
               send_mesg.push_back(d_local_max_box_size[bn][i]);
            }
         }
         send_mesg.push_back(getLocalNumberOfBoxes());
         send_mesg.push_back(-static_cast<int>(getLocalNumberOfBoxes()));
         send_mesg.push_back(getLocalNumberOfCells());
         send_mesg.push_back(-static_cast<int>(getLocalNumberOfCells()));

         std::vector<int> recv_mesg(send_mesg.size());
         d_mpi.Allreduce(
            &send_mesg[0],
            &recv_mesg[0],
            static_cast<int>(send_mesg.size()),
            MPI_INT,
            MPI_MAX);

         int tmpi = -1;
         for (int bn = 0; bn < nblocks; ++bn) {
            for (int i = 0; i < dim.getValue(); ++i) {
               d_global_bounding_box[bn].lower()[i] = -recv_mesg[++tmpi];
               d_global_bounding_box[bn].upper()[i] = recv_mesg[++tmpi];
               d_global_min_box_size[bn][i] = -recv_mesg[++tmpi];
               d_global_max_box_size[bn][i] = recv_mesg[++tmpi];
            }
         }
         d_max_number_of_boxes = recv_mesg[++tmpi];
         d_min_number_of_boxes = -recv_mesg[++tmpi];
         d_max_number_of_cells = recv_mesg[++tmpi];
         d_min_number_of_cells = -recv_mesg[++tmpi];
         TBOX_ASSERT(tmpi == int(recv_mesg.size() - 1));

      } else {
         d_global_bounding_box = d_local_bounding_box;
      }
   }

   d_global_data_up_to_date = true;

   t_cache_global_reduced_data->stop();
}

int BoxLevel::getGlobalNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_boxes;
}

int BoxLevel::getMaxNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_max_number_of_boxes;
}

int BoxLevel::getMinNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_min_number_of_boxes;
}

int BoxLevel::getGlobalNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_cells;
}

int BoxLevel::getMaxNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_max_number_of_cells;
}

int BoxLevel::getMinNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_min_number_of_cells;
}

const Box& BoxLevel::getGlobalBoundingBox(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_bounding_box[block_num];
}

const Box& BoxLevel::getLocalBoundingBox(int block_num) const
{
   return d_local_bounding_box[block_num];
}

const IntVector& BoxLevel::getLocalMaxBoxSize(int block_num) const
{
   return d_local_max_box_size[block_num];
}

const IntVector& BoxLevel::getGlobalMaxBoxSize(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_max_box_size[block_num];
}

const IntVector& BoxLevel::getLocalMinBoxSize(int block_num) const
{
   return d_local_min_box_size[block_num];
}

const IntVector& BoxLevel::getGlobalMinBoxSize(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_min_box_size[block_num];
}

bool BoxLevel::getSpatiallyEqualBox(
   const Box& box_to_match,
   const BlockId& block_id,
   Box& matching_box) const
{
   bool box_exists = false;
   for (BoxContainerSingleBlockIterator itr(d_boxes, block_id);
        itr.isValid(); ++itr) {
      if (box_to_match.isSpatiallyEqual(*itr)) {
         box_exists = true;
         matching_box = *itr;
         break;
      }
   }
   return box_exists;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevel::setParallelState(
   const ParallelState parallel_state)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isInitialized()) {
      TBOX_ERROR(
         "BoxLevel::setParallelState: Cannot change the parallel state of\n"
         << "an uninitialized BoxLevel.  See BoxLevel::initialize()");
   }
#endif
   if (parallel_state != DISTRIBUTED && parallel_state != GLOBALIZED) {
      TBOX_ERROR(
         "BoxLevel::setParallelState: Invalid distribution state: "
         << parallel_state << "\n");
   }

   if (d_parallel_state == DISTRIBUTED && parallel_state == GLOBALIZED) {
      acquireRemoteBoxes();
   } else if (d_parallel_state == GLOBALIZED && parallel_state ==
              DISTRIBUTED) {
      d_global_boxes.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevel::acquireRemoteBoxes()
{
   BoxLevel* object = this;
   acquireRemoteBoxes(1, &object);
}

/*
 ***********************************************************************
 * Acquire remote Boxes for multiple BoxLevels.
 * This method combines communication for the multiple
 * BoxLevels to increase message passing efficiency.
 *
 * Note: This method is stateless (could be static).
 ***********************************************************************
 */

void BoxLevel::acquireRemoteBoxes(
   const int num_sets,
   BoxLevel* multiple_box_levels[])
{
   if (d_mpi.getSize() == 1) {
      // In single-proc mode, we already have all the Boxes already.
      for (int n = 0; n < num_sets; ++n) {
         multiple_box_levels[n]->d_global_boxes =
            multiple_box_levels[n]->d_boxes;
      }
      return;
   }

   t_acquire_remote_boxes->start();
   int n;

#ifdef DEBUG_CHECK_ASSERTIONS
   for (n = 0; n < num_sets; ++n) {
      if (multiple_box_levels[n]->getParallelState() !=
          DISTRIBUTED) {
         TBOX_ERROR("BoxLevel objects must be in distributed mode\n"
            << "when acquiring remote boxes.\n");
      }
   }
#endif

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack Boxes from all BoxLevels into a single message.
    */
   for (n = 0; n < num_sets; ++n) {
      const BoxLevel& box_level =
         *multiple_box_levels[n];
      box_level.acquireRemoteBoxes_pack(send_mesg);
   }
   int send_mesg_size = static_cast<int>(send_mesg.size());

   /*
    * Send and receive the data.
    */

   std::vector<int> recv_mesg_size(d_mpi.getSize());
   d_mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(d_mpi.getSize());
   int totl_size = 0;
   for (n = 0; n < d_mpi.getSize(); ++n) {
      proc_offset[n] = totl_size;
      totl_size += recv_mesg_size[n];
   }
   recv_mesg.resize(totl_size, BAD_INT);
   d_mpi.Allgatherv(&send_mesg[0],
      send_mesg_size,
      MPI_INT,
      &recv_mesg[0],
      &recv_mesg_size[0],
      &proc_offset[0],
      MPI_INT);

   /*
    * Extract Box info received from other processors.
    */
   for (n = 0; n < num_sets; ++n) {
      BoxLevel& box_level =
         *multiple_box_levels[n];
      box_level.acquireRemoteBoxes_unpack(recv_mesg,
         proc_offset);
   }

   t_acquire_remote_boxes->stop();

}

/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevel::acquireRemoteBoxes_pack(
   std::vector<int>& send_mesg) const
{
   const tbox::Dimension& dim(getDim());
   /*
    * Box acquisition occurs during globalization.  Thus, do not
    * rely on current value of d_parallel_state.
    */

   /*
    * Pack Box info from d_boxes into send_mesg,
    * starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of Boxes from self
    *   - Self Boxes
    */
   const int box_com_buf_size = Box::commBufferSize(dim);
   const int send_mesg_size = 1 + box_com_buf_size
      * static_cast<int>(d_boxes.size());
   const int old_size = static_cast<int>(send_mesg.size());
   send_mesg.resize(old_size + send_mesg_size, BAD_INT);

   int* ptr = &send_mesg[0] + old_size;
   *(ptr++) = static_cast<int>(d_boxes.size());

   for (BoxContainer::ConstIterator i_boxes = d_boxes.begin();
        i_boxes != d_boxes.end();
        ++i_boxes) {
      (*i_boxes).putToIntBuffer(ptr);
      ptr += box_com_buf_size;
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevel::acquireRemoteBoxes_unpack(
   const std::vector<int>& recv_mesg,
   std::vector<int>& proc_offset)
{
   const tbox::Dimension& dim(getDim());
   /*
    * Unpack Box info from recv_mesg into d_global_boxes,
    * starting at the offset location.
    * Advance the proc_offset past the used data.
    */
   int n;
   int box_com_buf_size = Box::commBufferSize(dim);

   for (n = 0; n < d_mpi.getSize(); ++n) {
      if (n != d_mpi.getRank()) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int n_self_boxes = *(ptr++);
         proc_offset[d_mpi.getRank()] += (n_self_boxes) * box_com_buf_size;

         int i;
         Box box(dim);

         for (i = 0; i < n_self_boxes; ++i) {
            box.getFromIntBuffer(ptr);
            d_global_boxes.insert(
               d_global_boxes.end(), box);
            ptr += box_com_buf_size;
         }

      } else {
         d_global_boxes.insert(
            d_boxes.begin(), d_boxes.end());
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

BoxContainer::ConstIterator BoxLevel::addBox(
   const Box& box,
   const BlockId& block_id,
   const bool use_vacant_index)
{
   const tbox::Dimension& dim(getDim());
   /*
    * FIXME: bug: if some procs add a Box and others do not,
    * their d_computed_global_* flags will be inconsistent
    * resulting in incomplete participation in future communication
    * calls to compute those parameters.
    *
    * This problem is not exclusive to box adding.
    * It would also happen if some call initialize()
    * and others do not.  But because box adding is finer grained
    * than initialize(), it is more likely that some processors
    * will skip over the box adding.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually adding Boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state.");
   }
#endif

   clearForBoxChanges(false);

   BoxContainer::Iterator new_iterator(d_boxes);

   if (d_boxes.size() == 0) {
      Box new_box(
         box,
         LocalId::getZero(),
         d_mpi.getRank(),
         block_id,
         PeriodicShiftCatalog::getCatalog(dim)->getZeroShiftNumber());
      new_iterator = d_boxes.insert(d_boxes.end(), new_box);
   } else {
      // Set new_index to one more than the largest index used.
      BoxContainer::Iterator ni = d_boxes.end();
      do {
         TBOX_ASSERT(ni != d_boxes.begin());   // There should not be all periodic images.
         --ni;
      } while (ni->isPeriodicImage());
      LocalId new_index = ni->getLocalId() + 1;
      if (use_vacant_index) {
         TBOX_ASSERT(new_index >= 0);

         if (new_index.getValue() !=
             static_cast<int>(d_local_number_of_boxes)) {
            /*
             * There is a smaller unused index we can use for the new index.
             */
            for (new_index = 0, ni = d_boxes.begin();
                 ni != d_boxes.end();
                 ++ni) {
               if ((*ni).getBlockId() == block_id) {
                  if (new_index != (*ni).getLocalId()) {
                     break;
                  }
                  if (!ni->isPeriodicImage()) {
                     ++new_index;
                  }
               }
            }
            // We should have found an unused index.
            TBOX_ASSERT(ni != d_boxes.end());
         }
      }

      const Box new_box(
         box, new_index, d_mpi.getRank(), block_id);
      new_iterator = d_boxes.insert(ni, new_box);
   }

   const IntVector box_size(box.numberCells());
   ++d_local_number_of_boxes;
   d_local_number_of_cells += box.size();
   d_local_bounding_box[block_id.getBlockValue()] += box;
   d_local_max_box_size[block_id.getBlockValue()].max(box_size);
   d_local_min_box_size[block_id.getBlockValue()].min(box_size);
   d_global_data_up_to_date = false;

   return new_iterator;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addPeriodicBox(
   const Box& ref_box,
   const PeriodicId& shift_number)
{
   // FIXME: We don't allow individually adding remote Boxes even in globalized state.  We probably shouldn't allow adding remote images either.
#ifdef DEBUG_CHECK_ASSERTIONS
   if (shift_number ==
       PeriodicShiftCatalog::getCatalog(getDim())->getZeroShiftNumber()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox cannot be used to add regular box.");
   }
   if (d_parallel_state != GLOBALIZED && ref_box.getOwnerRank() !=
       d_mpi.getRank()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox: Cannot add remote Box\n"
         << "(owned by rank " << ref_box.getOwnerRank() << ")\n"
         << "when not in GLOBALIZED state.");
   }
#endif

   clearForBoxChanges(false);

   Box image_box(ref_box, shift_number, d_ratio);

#ifdef DEBUG_CHECK_ASSERTIONS
   BoxContainer& boxes =
      d_parallel_state == DISTRIBUTED ? d_boxes : d_global_boxes;
   /*
    * Sanity checks:
    *
    * - Require that the real version of the reference Box exists
    *   before adding the periodic image Box.
    */
   Box real_box(getDim(),
                       ref_box.getLocalId(),
                       ref_box.getOwnerRank(),
                       ref_box.getBlockId(),
                       PeriodicShiftCatalog::getCatalog(
                          getDim())->getZeroShiftNumber());
   if (boxes.find(real_box) == boxes.end()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox: cannot add periodic image Box "
         << image_box
         << "\nwithout the real Box (" << real_box
         << ") already in the BoxLevel.\n");
   }
#endif

   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(image_box);
   }
   if (image_box.getOwnerRank() == d_mpi.getRank()) {
      d_boxes.insert(image_box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addBox(
   const Box& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != GLOBALIZED && box.getOwnerRank() != d_mpi.getRank()) {
      TBOX_ERROR("BoxLevel::addBox: Cannot add remote Box\n"
         << "(owned by rank " << box.getOwnerRank() << ")\n"
         << "when not in GLOBALIZED state.");
   }
#endif

   clearForBoxChanges(false);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Sanity checks:
    * - Require that the real Box exists before adding the periodic image Box.
    */
   if (box.isPeriodicImage()) {
      Box real_box(getDim(),
                          box.getLocalId(),
                          box.getOwnerRank(),
                          box.getBlockId(),
                          PeriodicShiftCatalog::getCatalog(
                             getDim())->getZeroShiftNumber());
      BoxContainer& boxes = box.getOwnerRank() ==
         d_mpi.getRank() ? d_boxes : d_global_boxes;
      if (boxes.find(real_box) == boxes.end()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add periodic image Box "
            << box
            << "\nwithout the real Box (" << real_box
            << ") already in the BoxLevel.\n");
      }
      if (d_global_boxes.find(box) !=
          d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add Box "
            << box
            << "\nbecause it already exists ("
            << *boxes.find(box) << "\n");
      }
   }
#endif

   // Update counters.
   if (!box.isPeriodicImage()) {
      if (box.getOwnerRank() == d_mpi.getRank()) {
         const IntVector box_size(box.numberCells());
         ++d_local_number_of_boxes;
         d_local_number_of_cells += box.size();
         d_local_bounding_box[box.getBlockId().getBlockValue()] += box;
         d_local_max_box_size[box.getBlockId().getBlockValue()].max(box_size);
         d_local_min_box_size[box.getBlockId().getBlockValue()].min(box_size);
      }
      d_global_data_up_to_date = false;
      /*
       * FIXME: bug: if some procs add a real Box and others do not,
       * their d_global_data_up_to_date flags will be inconsistent
       * resulting in incomplete participation in future collective
       * communication to compute that parameter.
       */
   }

   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(box);
   }
   if (box.getOwnerRank() == d_mpi.getRank()) {
      d_boxes.insert(box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addBoxWithoutUpdate(
   const Box& box)
{
   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(box);
   }
   d_boxes.insert(box);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   BoxContainer::Iterator& ibox)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually erasing boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state.");
   }
#endif

   clearForBoxChanges();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (ibox != d_boxes.find(*ibox)) {
      TBOX_ERROR("BoxLevel::eraseBox: Attempt to erase a\n"
         << "Box that does not belong to the BoxLevel\n"
         << "object.\n");
   }
#endif

   if (ibox->isPeriodicImage()) {
      d_boxes.erase(ibox++);
      // No need to update counters (they neglect periodic images).
   } else {
      /*
       * Update counters.  Bounding box cannot be updated (without
       * recomputing) because we don't know how the erased Box
       * affects the bounding box.
       */
      d_local_bounding_box_up_to_date = d_global_data_up_to_date = false;
      --d_local_number_of_boxes;
      d_local_number_of_cells -= ibox->size();
      // Erase real Box and its periodic images.
      const LocalId& local_id = ibox->getLocalId();
      do {
         d_boxes.erase(ibox++);
      } while (ibox != d_boxes.end() && ibox->getLocalId() ==
               local_id);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   const Box& box)
{
   /*
    * FIXME: bug: if some procs erase some Boxes and others do
    * not, their d_computed_global_* flags will be inconsistent
    * resulting in incomplete participation in future communication
    * calls to compute those parameters.
    *
    * This problem is not exclusive to box adding/erasing.  It would
    * also happen if some call initialize() and others do not.  But
    * because box adding is finer grained than initialize(), it is
    * more likely that some processors will skip over the box adding.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually erasing Boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state.");
   }
#endif

   clearForBoxChanges();

   d_local_bounding_box_up_to_date = d_global_data_up_to_date = false;

   BoxContainer::Iterator ibox = d_boxes.find(box);
   if (ibox == d_boxes.end()) {
      TBOX_ERROR("BoxLevel::eraseBox: Box to be erased ("
         << box << ") is NOT a part of the BoxLevel.\n");
   }
   d_boxes.erase(ibox);
}

/*
 ****************************************************************************
 ****************************************************************************
 */
void BoxLevel::eraseBoxWithoutUpdate(
   const Box& box)
{
   d_boxes.erase(box);
   return;
}

/*
 ****************************************************************************
 ****************************************************************************
 */
void BoxLevel::refineBoxes(
   BoxLevel& finer,
   const IntVector& ratio,
   const IntVector& final_ratio) const
{
   finer.d_boxes = d_boxes;
   finer.d_boxes.refine(ratio);
   finer.d_ratio = final_ratio;
   return;
}

/*
 ****************************************************************************
 ****************************************************************************
 */
void BoxLevel::coarsenBoxes(
   BoxLevel& coarser,
   const IntVector& ratio,
   const IntVector& final_ratio) const
{
   coarser.d_boxes = d_boxes;
   coarser.d_boxes.coarsen(ratio);
   coarser.d_ratio = final_ratio;
   return;
}

/*
 ****************************************************************************
 ****************************************************************************
 */
const BoxLevel& BoxLevel::getGlobalizedVersion() const
{
   TBOX_ASSERT(isInitialized());

   if (d_parallel_state == GLOBALIZED) {
      return *this;
   }

   if (d_globalized_version == NULL) {
      BoxLevel* globalized_version = new BoxLevel(*this);
      globalized_version->setParallelState(GLOBALIZED);
      TBOX_ASSERT(globalized_version->getParallelState() == GLOBALIZED);
      d_globalized_version = globalized_version;
      globalized_version = NULL;
   }

   TBOX_ASSERT(d_globalized_version->getParallelState() == GLOBALIZED);
   return *d_globalized_version;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
PersistentOverlapConnectors& BoxLevel::getPersistentOverlapConnectors()
const
{
   if (d_persistent_overlap_connectors == NULL) {
      d_persistent_overlap_connectors = new PersistentOverlapConnectors(*this);
   }
   return *d_persistent_overlap_connectors;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void BoxLevel::getGlobalBoxes(BoxContainer& global_boxes) const
{
   for (BoxContainer::ConstIterator itr = d_global_boxes.begin();
        itr != d_global_boxes.end(); ++itr) {
      global_boxes.pushBack(*itr);
   }
}

/*
 ***********************************************************************
 * Write the BoxLevel to a database.
 *
 * Write only local parts.
 ***********************************************************************
 */

void BoxLevel::putToDatabase(
   tbox::Database& database) const
{
   database.putBool("d_is_mapped_box_level", true);
   database.putInteger(
      "HIER_MAPPED_BOX_LEVEL_VERSION", HIER_MAPPED_BOX_LEVEL_VERSION);
   database.putInteger("d_nproc", d_mpi.getSize());
   database.putInteger("d_rank", d_mpi.getRank());
   database.putInteger("dim", d_ratio.getDim().getValue());
   database.putIntegerArray("d_ratio", &d_ratio[0], d_ratio.getDim().getValue());
   getBoxes().putToDatabase(*database.putDatabase("mapped_boxes"));
}

/*
 ***********************************************************************
 * Read the BoxLevel from a database.
 ***********************************************************************
 */

void BoxLevel::getFromDatabase(
   tbox::Database& database,
   const tbox::ConstPointer<GridGeometry>& grid_geom)
{
   TBOX_ASSERT(database.isInteger("dim"));
   const tbox::Dimension dim(static_cast<unsigned short>(database.getInteger("dim")));
   TBOX_ASSERT(getDim() == dim);

   IntVector ratio(dim);
   database.getIntegerArray("d_ratio", &ratio[0], dim.getValue());

#ifdef DEBUG_CHECK_ASSERTIONS
   const int version = database.getInteger("HIER_MAPPED_BOX_LEVEL_VERSION");
   const int nproc = database.getInteger("d_nproc");
   const int rank = database.getInteger("d_rank");
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));
   TBOX_ASSERT(version <= HIER_MAPPED_BOX_LEVEL_VERSION);
#endif

   /*
    * If the communicator is already set, use it.  Otherwise, use the
    * one in tbox::SAMRAI_MPI::getSAMRAIWorld().
    *
    * There must be a better way to handle the communicator than this.
    */
   if (isInitialized()) {
      TBOX_ASSERT(ratio == d_ratio);
      if (d_parallel_state != DISTRIBUTED) {
         setParallelState(DISTRIBUTED);
         d_boxes.clear();
      }
   } else {
      TBOX_WARNING(
         "BoxLevel::getFromDatabase: Uninitialized MPI communicator.\n"
         << "Using tbox::SAMRAI_MPI::getSAMRAIWorld().");
      initialize(ratio, grid_geom,
         tbox::SAMRAI_MPI::getSAMRAIWorld(), DISTRIBUTED);
   }

   /*
    * Failing these asserts means that we don't have a compatible
    * database for the number of processors or we are reading another
    * processor's data.
    */
   TBOX_ASSERT(nproc == d_mpi.getSize());
   TBOX_ASSERT(rank == d_mpi.getRank());

   d_boxes.getFromDatabase(*database.getDatabase("mapped_boxes"));
   computeLocalRedundantData();

}

/*
 ***********************************************************************
 * Construct a BoxLevel Outputter with formatting parameters.
 ***********************************************************************
 */

BoxLevel::Outputter::Outputter(
   const BoxLevel& box_level,
   const std::string& border,
   int detail_depth):
   d_level(box_level),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a BoxLevel according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
   std::ostream& s,
   const BoxLevel::Outputter& format)
{
   format.d_level.recursivePrint(s, format.d_border, format.d_detail_depth);
   return s;
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxLevel to a stream.
 ***********************************************************************
 */

BoxLevel::Outputter BoxLevel::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}

/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 * Print out global bounding box only if it has been computed already.
 ***********************************************************************
 */

void BoxLevel::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   if (detail_depth < 0) return;

   if (!isInitialized()) {
      co << border << "Uninitialized.\n";
      return;
   }
   co // << "Address        : " << (void*)this << '\n'
   << border << "Parallel state : "
   << (getParallelState() == DISTRIBUTED ? "DIST" : "GLOB") << '\n'
   << border << "Ratio          : " << getRefinementRatio() << '\n'
   << border << "Box count      : " << d_local_number_of_boxes << ", "
   << d_global_number_of_boxes << '\n'
   << border << "Cell count     : " << d_local_number_of_cells << ", "
   << d_global_number_of_cells << '\n'
   << border << "Bounding box   : " << getLocalBoundingBox(0) << ", "
   << (d_global_data_up_to_date ? getGlobalBoundingBox(0) : Box(getDim()))
   << '\n'
   << border << "Comm,rank,nproc: " << d_mpi.getCommunicator() << ", " << d_mpi.getRank()
   << ", " << d_mpi.getSize() << '\n'
   ;
   if (detail_depth > 0) {
      co << border << "Boxes:\n";
      if (getParallelState() == GLOBALIZED) {
         /*
          * Print boxes from all ranks.
          */
         for (BoxContainer::ConstIterator bi = d_global_boxes.begin();
              bi != d_global_boxes.end();
              ++bi) {
            Box box = *bi;
            co << border << "    "
            << box << "   "
            << box.numberCells() << '\n';
         }
      } else {
         /*
          * Print local boxes only.
          */
         for (BoxContainer::ConstIterator bi = d_boxes.begin();
              bi != d_boxes.end();
              ++bi) {
            Box box = *bi;
            co << border << "    "
            << box << "   "
            << box.numberCells() << '\n';
         }
      }
   }
}


/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevel::initializeCallback()
{
   t_initialize_private = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::initializePrivate()");
   t_acquire_remote_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::acquireRemoteBoxes()");
   t_cache_global_reduced_data = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::cacheGlobalReducedData()");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */

void BoxLevel::finalizeCallback()
{
   t_initialize_private.setNull();
   t_acquire_remote_boxes.setNull();
   t_cache_global_reduced_data.setNull();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
