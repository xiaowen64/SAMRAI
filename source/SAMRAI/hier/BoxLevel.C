/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of mapped_boxes in a mapped_box_level of a distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevel_C
#define included_hier_BoxLevel_C

#include "SAMRAI/hier/BoxLevel.h"

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
tbox::Pointer<tbox::Timer> BoxLevel::t_acquire_remote_mapped_boxes;
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
   d_mapped_boxes(tbox::Dimension::getInvalidDimension()),
   d_global_mapped_boxes(tbox::Dimension::getInvalidDimension()),
   d_ratio(tbox::Dimension::getInvalidDimension(), 0),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_mapped_boxes(0),
   d_global_number_of_mapped_boxes(-1),

   d_max_number_of_mapped_boxes(-1),
   d_min_number_of_mapped_boxes(-1),
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
   d_mapped_boxes(dim),
   d_global_mapped_boxes(dim),
   d_ratio(dim, 0),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_mapped_boxes(0),
   d_global_number_of_mapped_boxes(-1),

   d_max_number_of_mapped_boxes(-1),
   d_min_number_of_mapped_boxes(-1),
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
   d_mapped_boxes(rhs.d_mapped_boxes),
   d_global_mapped_boxes(rhs.d_global_mapped_boxes),
   d_ratio(rhs.d_ratio),

   d_local_number_of_cells(rhs.d_local_number_of_cells),
   d_global_number_of_cells(rhs.d_global_number_of_cells),
   d_local_number_of_mapped_boxes(rhs.d_local_number_of_mapped_boxes),
   d_global_number_of_mapped_boxes(rhs.d_global_number_of_mapped_boxes),

   d_max_number_of_mapped_boxes(rhs.d_max_number_of_mapped_boxes),
   d_min_number_of_mapped_boxes(rhs.d_min_number_of_mapped_boxes),
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
   const BoxSet& mapped_boxes,
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state):
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mapped_boxes(ratio.getDim()),
   d_global_mapped_boxes(ratio.getDim()),
   d_ratio(ratio),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_mapped_boxes(0),
   d_global_number_of_mapped_boxes(-1),

   d_max_number_of_mapped_boxes(-1),
   d_min_number_of_mapped_boxes(-1),
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
   initialize(mapped_boxes, ratio, grid_geom, mpi, parallel_state);
}

BoxLevel::BoxLevel(
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state):
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mapped_boxes(ratio.getDim()),
   d_global_mapped_boxes(ratio.getDim()),
   d_ratio(ratio),

   d_local_number_of_cells(0),
   d_global_number_of_cells(-1),
   d_local_number_of_mapped_boxes(0),
   d_global_number_of_mapped_boxes(-1),

   d_max_number_of_mapped_boxes(-1),
   d_min_number_of_mapped_boxes(-1),
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
   BoxSet dummy_mapped_boxes(ratio.getDim());
   initialize(dummy_mapped_boxes, ratio, grid_geom, mpi, parallel_state);
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
   const BoxSet& mapped_boxes,
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   if (&mapped_boxes != &d_mapped_boxes) {
      d_mapped_boxes = mapped_boxes;
   }
   initializePrivate(ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void BoxLevel::initialize(
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   d_mapped_boxes.clear();
   initializePrivate(
      ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void BoxLevel::swapInitialize(
   BoxSet& mapped_boxes,
   const IntVector& ratio,
   const tbox::ConstPointer<GridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   TBOX_ASSERT(&mapped_boxes != &d_mapped_boxes);   // Library error if this fails.
   d_mapped_boxes.swap(mapped_boxes);
   initializePrivate(ratio,
      grid_geom,
      mpi,
      parallel_state);
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
      d_global_mapped_boxes.clear();
   } else {
      d_global_mapped_boxes = d_mapped_boxes;
   }

   // Erase non-local Boxes, if any, from d_mapped_boxes.
   for (BoxSet::SetIterator mbi(d_mapped_boxes.setBegin());
        mbi != d_mapped_boxes.setEnd(); /* incremented in loop */) {
      if (mbi->getOwnerRank() != d_mpi.getRank()) {
         d_mapped_boxes.erase(mbi++);
      } else {
         ++mbi;
      }
   }

   d_ratio = ratio;
   d_parallel_state = parallel_state;
   d_global_number_of_cells = -1;
   d_global_number_of_mapped_boxes = -1;
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
      d_mapped_boxes.removePeriodicImageBoxes();
      if ( d_parallel_state == GLOBALIZED ) {
         d_global_mapped_boxes.removePeriodicImageBoxes();
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
      d_mapped_boxes.clear();
      d_global_mapped_boxes.clear();
      d_ratio(0) = 0;
      d_local_number_of_cells = 0;
      d_global_number_of_cells = -1;
      d_local_number_of_mapped_boxes = 0;
      d_global_number_of_mapped_boxes = -1;
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
      level_a.d_mapped_boxes.swap(level_b.d_mapped_boxes);
      level_a.d_global_mapped_boxes.swap(level_b.d_global_mapped_boxes);
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

      tmpint = static_cast<int>(level_a.d_local_number_of_mapped_boxes);
      level_a.d_local_number_of_mapped_boxes = level_b.d_local_number_of_mapped_boxes;
      level_b.d_local_number_of_mapped_boxes = tmpint;

      tmpint = level_a.d_global_number_of_mapped_boxes;
      level_a.d_global_number_of_mapped_boxes = level_b.d_global_number_of_mapped_boxes;
      level_b.d_global_number_of_mapped_boxes = tmpint;

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

   d_local_number_of_mapped_boxes = 0;
   d_local_number_of_cells = 0;

   if (int(d_local_bounding_box.size()) != nblocks) {
      d_local_bounding_box.resize(nblocks, hier::Box(d_grid_geometry->getDim()));
      d_local_min_box_size.resize(nblocks, max_vec);
      d_local_max_box_size.resize(nblocks, zero_vec);
   }

   for (RealBoxConstIterator ni(d_mapped_boxes); ni.isValid(); ++ni) {

      int block_num = ni->getBlockId().getBlockValue();
      const IntVector boxdim(ni->numberCells());
      ++d_local_number_of_mapped_boxes;
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
      d_global_number_of_mapped_boxes = 0;
      d_global_number_of_cells = 0;
      for (RealBoxConstIterator ni(d_global_mapped_boxes);
           ni.isValid();
           ++ni) {
         ++d_global_number_of_mapped_boxes;
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
         d_global_number_of_mapped_boxes = tmpb[0];
         d_global_number_of_cells = tmpb[1];
      } else {
         d_global_number_of_mapped_boxes = getLocalNumberOfBoxes();
         d_global_number_of_cells = getLocalNumberOfCells();
      }

      TBOX_ASSERT(d_global_number_of_mapped_boxes >= 0);
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
      d_max_number_of_mapped_boxes = d_min_number_of_mapped_boxes = getLocalNumberOfBoxes();
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
         d_max_number_of_mapped_boxes = recv_mesg[++tmpi];
         d_min_number_of_mapped_boxes = -recv_mesg[++tmpi];
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
   return d_global_number_of_mapped_boxes;
}

int BoxLevel::getMaxNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_max_number_of_mapped_boxes;
}

int BoxLevel::getMinNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_min_number_of_mapped_boxes;
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
      d_global_mapped_boxes.clear();
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
   BoxLevel* multiple_mapped_box_levels[])
{
   if (d_mpi.getSize() == 1) {
      // In single-proc mode, we already have all the Boxes already.
      for (int n = 0; n < num_sets; ++n) {
         multiple_mapped_box_levels[n]->d_global_mapped_boxes =
            multiple_mapped_box_levels[n]->d_mapped_boxes;
      }
      return;
   }

   t_acquire_remote_mapped_boxes->start();
   int n;

#ifdef DEBUG_CHECK_ASSERTIONS
   for (n = 0; n < num_sets; ++n) {
      if (multiple_mapped_box_levels[n]->getParallelState() !=
          DISTRIBUTED) {
         TBOX_ERROR("BoxLevel objects must be in distributed mode\n"
            << "when acquiring remote mapped_boxes.\n");
      }
   }
#endif

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack Boxes from all BoxLevels into a single message.
    */
   for (n = 0; n < num_sets; ++n) {
      const BoxLevel& mapped_box_level =
         *multiple_mapped_box_levels[n];
      mapped_box_level.acquireRemoteBoxes_pack(send_mesg);
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
      BoxLevel& mapped_box_level =
         *multiple_mapped_box_levels[n];
      mapped_box_level.acquireRemoteBoxes_unpack(recv_mesg,
         proc_offset);
   }

   t_acquire_remote_mapped_boxes->stop();

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
    * Pack Box info from d_mapped_boxes into send_mesg,
    * starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of Boxes from self
    *   - Self Boxes
    */
   const int mapped_box_com_buf_size = Box::commBufferSize(dim);
   const int send_mesg_size = 1 + mapped_box_com_buf_size
      * static_cast<int>(d_mapped_boxes.size());
   const int old_size = static_cast<int>(send_mesg.size());
   send_mesg.resize(old_size + send_mesg_size, BAD_INT);

   int* ptr = &send_mesg[0] + old_size;
   *(ptr++) = static_cast<int>(d_mapped_boxes.size());

   for (BoxSet::SetConstIterator i_mapped_boxes = d_mapped_boxes.setBegin();
        i_mapped_boxes != d_mapped_boxes.setEnd();
        ++i_mapped_boxes) {
      (*i_mapped_boxes).putToIntBuffer(ptr);
      ptr += mapped_box_com_buf_size;
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
    * Unpack Box info from recv_mesg into d_global_mapped_boxes,
    * starting at the offset location.
    * Advance the proc_offset past the used data.
    */
   int n;
   int mapped_box_com_buf_size = Box::commBufferSize(dim);

   for (n = 0; n < d_mpi.getSize(); ++n) {
      if (n != d_mpi.getRank()) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int n_self_mapped_boxes = *(ptr++);
         proc_offset[d_mpi.getRank()] += (n_self_mapped_boxes) * mapped_box_com_buf_size;

         int i;
         Box mapped_box(dim);

         for (i = 0; i < n_self_mapped_boxes; ++i) {
            mapped_box.getFromIntBuffer(ptr);
            d_global_mapped_boxes.insert(
               d_global_mapped_boxes.setEnd(), mapped_box);
            ptr += mapped_box_com_buf_size;
         }

      } else {
         for (BoxContainer::SetConstIterator ni = d_mapped_boxes.setBegin();
              ni != d_mapped_boxes.setEnd(); ++ni) {
            d_global_mapped_boxes.insert(*ni);
         }
//         d_global_mapped_boxes.insert(
//            d_mapped_boxes.setBegin(), d_mapped_boxes.setEnd());
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

BoxSet::SetIterator BoxLevel::addBox(
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

   BoxSet::SetIterator new_iterator(d_mapped_boxes);

   if (d_mapped_boxes.size() == 0) {
      Box new_mapped_box =
         Box(box,
            LocalId::getZero(),
            d_mpi.getRank(),
            block_id,
            PeriodicShiftCatalog::getCatalog(dim)->getZeroShiftNumber());
      new_iterator = d_mapped_boxes.insert(d_mapped_boxes.setEnd(), new_mapped_box);
   } else {
      // Set new_index to one more than the largest index used.
      BoxSet::SetIterator ni = d_mapped_boxes.setEnd();
      do {
         TBOX_ASSERT(ni != d_mapped_boxes.setBegin());   // There should not be all periodic images.
         --ni;
      } while (ni->isPeriodicImage());
      LocalId new_index = ni->getLocalId() + 1;
      if (use_vacant_index) {
         TBOX_ASSERT(new_index >= 0);

         if (new_index.getValue() !=
             static_cast<int>(d_local_number_of_mapped_boxes)) {
            /*
             * There is a smaller unused index we can use for the new index.
             */
            for (new_index = 0, ni = d_mapped_boxes.setBegin();
                 ni != d_mapped_boxes.setEnd();
                 ++ni) {
               if (new_index != (*ni).getLocalId()) {
                  break;
               }
               if (!ni->isPeriodicImage()) {
                  ++new_index;
               }
            }
            // We should have found an unused index.
            TBOX_ASSERT(ni != d_mapped_boxes.setEnd());
         }
      }

      const Box new_mapped_box =
         Box(box, new_index, d_mpi.getRank(), block_id);
      new_iterator = d_mapped_boxes.insert(ni, new_mapped_box);
   }

   const IntVector box_size(box.numberCells());
   ++d_local_number_of_mapped_boxes;
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
   const Box& ref_mapped_box,
   const PeriodicId& shift_number)
{
   // FIXME: We don't allow individually adding remote Boxes even in globalized state.  We probably shouldn't allow adding remote images either.
#ifdef DEBUG_CHECK_ASSERTIONS
   if (shift_number ==
       PeriodicShiftCatalog::getCatalog(getDim())->getZeroShiftNumber()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox cannot be used to add regular mapped_box.");
   }
   if (d_parallel_state != GLOBALIZED && ref_mapped_box.getOwnerRank() !=
       d_mpi.getRank()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox: Cannot add remote Box\n"
         << "(owned by rank " << ref_mapped_box.getOwnerRank() << ")\n"
         << "when not in GLOBALIZED state.");
   }
#endif

   clearForBoxChanges(false);

   Box image_mapped_box(ref_mapped_box, shift_number, d_ratio);

#ifdef DEBUG_CHECK_ASSERTIONS
   BoxSet& mapped_boxes =
      d_parallel_state == DISTRIBUTED ? d_mapped_boxes : d_global_mapped_boxes;
   /*
    * Sanity checks:
    *
    * - Require that the real version of the reference Box exists
    *   before adding the periodic image Box.
    */
   Box real_mapped_box(getDim(),
                       ref_mapped_box.getLocalId(),
                       ref_mapped_box.getOwnerRank(),
                       ref_mapped_box.getBlockId(),
                       PeriodicShiftCatalog::getCatalog(
                          getDim())->getZeroShiftNumber());
   if (mapped_boxes.find(real_mapped_box) == mapped_boxes.setEnd()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox: cannot add periodic image Box "
         << image_mapped_box
         << "\nwithout the real Box (" << real_mapped_box
         << ") already in the BoxLevel.\n");
   }
#endif

   if (d_parallel_state == GLOBALIZED) {
      d_global_mapped_boxes.insert(image_mapped_box);
   }
   if (image_mapped_box.getOwnerRank() == d_mpi.getRank()) {
      d_mapped_boxes.insert(image_mapped_box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addBox(
   const Box& mapped_box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != GLOBALIZED && mapped_box.getOwnerRank() != d_mpi.getRank()) {
      TBOX_ERROR("BoxLevel::addBox: Cannot add remote Box\n"
         << "(owned by rank " << mapped_box.getOwnerRank() << ")\n"
         << "when not in GLOBALIZED state.");
   }
#endif

   clearForBoxChanges(false);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Sanity checks:
    * - Require that the real Box exists before adding the periodic image Box.
    */
   if (mapped_box.isPeriodicImage()) {
      Box real_mapped_box(getDim(),
                          mapped_box.getLocalId(),
                          mapped_box.getOwnerRank(),
                          mapped_box.getBlockId(),
                          PeriodicShiftCatalog::getCatalog(
                             getDim())->getZeroShiftNumber());
      BoxSet& mapped_boxes = mapped_box.getOwnerRank() ==
         d_mpi.getRank() ? d_mapped_boxes : d_global_mapped_boxes;
      if (mapped_boxes.find(real_mapped_box) == mapped_boxes.setEnd()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add periodic image Box "
            << mapped_box
            << "\nwithout the real Box (" << real_mapped_box
            << ") already in the BoxLevel.\n");
      }
      if (d_global_mapped_boxes.find(mapped_box) !=
          d_global_mapped_boxes.setEnd()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add Box "
            << mapped_box
            << "\nbecause it already exists ("
            << *mapped_boxes.find(mapped_box) << "\n");
      }
   }
#endif

   // Update counters.
   if (!mapped_box.isPeriodicImage()) {
      if (mapped_box.getOwnerRank() == d_mpi.getRank()) {
         const IntVector box_size(mapped_box.numberCells());
         ++d_local_number_of_mapped_boxes;
         d_local_number_of_cells += mapped_box.size();
         d_local_bounding_box[mapped_box.getBlockId().getBlockValue()] += mapped_box;
         d_local_max_box_size[mapped_box.getBlockId().getBlockValue()].max(box_size);
         d_local_min_box_size[mapped_box.getBlockId().getBlockValue()].min(box_size);
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
      d_global_mapped_boxes.insert(mapped_box);
   }
   if (mapped_box.getOwnerRank() == d_mpi.getRank()) {
      d_mapped_boxes.insert(mapped_box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   BoxSet::SetIterator& ibox)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually erasing mapped_boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state.");
   }
#endif

   clearForBoxChanges();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (ibox != d_mapped_boxes.find(*ibox)) {
      TBOX_ERROR("BoxLevel::eraseBox: Attempt to erase a\n"
         << "Box that does not belong to the BoxLevel\n"
         << "object.\n");
   }
#endif

   if (ibox->isPeriodicImage()) {
      d_mapped_boxes.erase(ibox++);
      // No need to update counters (they neglect periodic images).
   } else {
      /*
       * Update counters.  Bounding box cannot be updated (without
       * recomputing) because we don't know how the erased Box
       * affects the bounding box.
       */
      d_local_bounding_box_up_to_date = d_global_data_up_to_date = false;
      --d_local_number_of_mapped_boxes;
      d_local_number_of_cells -= ibox->size();
      // Erase real Box and its periodic images.
      const LocalId& local_id = ibox->getLocalId();
      do {
         d_mapped_boxes.erase(ibox++);
      } while (ibox != d_mapped_boxes.setEnd() && ibox->getLocalId() ==
               local_id);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   const Box& mapped_box)
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

   BoxSet::SetIterator ibox = d_mapped_boxes.find(mapped_box);
   if (ibox == d_mapped_boxes.setEnd()) {
      TBOX_ERROR("BoxLevel::eraseBox: Box to be erased ("
         << mapped_box << ") is NOT a part of the BoxLevel.\n");
   }
   d_mapped_boxes.erase(ibox);
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
void BoxLevel::getGlobalBoxes(BoxList& global_boxes) const
{
   for (BoxSet::SetConstIterator itr = d_global_mapped_boxes.setBegin();
        itr != d_global_mapped_boxes.setEnd(); itr++) {
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
         d_mapped_boxes.clear();
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

   d_mapped_boxes.getFromDatabase(*database.getDatabase("mapped_boxes"));
   computeLocalRedundantData();

}

/*
 ***********************************************************************
 * Construct a BoxLevel Outputter with formatting parameters.
 ***********************************************************************
 */

BoxLevel::Outputter::Outputter(
   const BoxLevel& mapped_box_level,
   const std::string& border,
   int detail_depth):
   d_level(mapped_box_level),
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
   << border << "Box count      : " << d_local_number_of_mapped_boxes << ", "
   << d_global_number_of_mapped_boxes << '\n'
   << border << "Cell count     : " << d_local_number_of_cells << ", "
   << d_global_number_of_cells << '\n'
   << border << "Bounding box   : " << getLocalBoundingBox(0) << ", "
   << (d_global_data_up_to_date ? getGlobalBoundingBox(0) : Box(getDim()))
   << '\n'
   << border << "Comm,rank,nproc: " << d_mpi.getCommunicator() << ", " << d_mpi.getRank()
   << ", " << d_mpi.getSize() << '\n'
   ;
   if (detail_depth > 0) {
      co << border << "Mapped_boxes:\n";
      if (getParallelState() == GLOBALIZED) {
         /*
          * Print mapped_boxes from all ranks.
          */
         for (BoxSet::SetConstIterator bi = d_global_mapped_boxes.setBegin();
              bi != d_global_mapped_boxes.setEnd();
              ++bi) {
            Box mapped_box = *bi;
            co << border << "    "
            << mapped_box << "   "
            << mapped_box.numberCells() << '\n';
         }
      } else {
         /*
          * Print local mapped_boxes only.
          */
         for (BoxSet::SetConstIterator bi = d_mapped_boxes.setBegin();
              bi != d_mapped_boxes.setEnd();
              ++bi) {
            Box mapped_box = *bi;
            co << border << "    "
            << mapped_box << "   "
            << mapped_box.numberCells() << '\n';
         }
      }
   }
}

/*
 ***********************************************************************
 * Write out some statistics on the mapped_boxes, including statistics
 * for judging mesh quality.
 ***********************************************************************
 */

void BoxLevel::printBoxStats(
   std::ostream& co,
   const std::string& border) const
{
   if (!isInitialized()) {
      co << "Uninitialized.\n";
      return;
   }

   const tbox::Dimension& dim(getDim());

   cacheGlobalReducedData();

   double ideal_surfarea =
      pow((double)getGlobalNumberOfCells() / d_mpi.getSize(), double(dim.getValue() - 1) / dim.getValue());

   // Per-processor statistics.
   double has_mapped_box = (double)(getLocalNumberOfBoxes() > 0);
   double number_of_mapped_boxes = (double)(getLocalNumberOfBoxes());
   double largest_dim = 0;
   double smallest_dim = (double)(getLocalNumberOfBoxes() == 0 ? 0 : 9999999);
   double largest_aspect = 0;
   double smallest_aspect = 1.;
   double sum_aspect = 0.;
   double sum_surfarea = 0.;
   double sum_normsurfarea = 0.;

   const BoxSet& mapped_boxes = getBoxes();

   for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {

      const Box& mapped_box = *ni;
      const IntVector boxdims = mapped_box.numberCells();
      const double boxvol = boxdims.getProduct();
      const double longdim = boxdims.max();
      const double shortdim = boxdims.min();
      const double aspect = (double)longdim / shortdim;
      double surfarea = 0.;
      for (int d = 0; d < dim.getValue(); ++d) {
         surfarea += 2 * (double)boxvol / boxdims(d);
      }

      largest_dim = tbox::MathUtilities<double>::Max(largest_dim, longdim);
      smallest_dim = tbox::MathUtilities<double>::Min(smallest_dim, shortdim);

      largest_aspect = tbox::MathUtilities<double>::Max(largest_aspect, aspect);
      smallest_aspect = tbox::MathUtilities<double>::Min(smallest_aspect,
            aspect);
      sum_aspect += aspect;

      sum_surfarea += surfarea;

   }

   sum_normsurfarea = sum_surfarea / ideal_surfarea;

   /*
    * Collect local statistics into arrays for global reduction.
    */
   {
      // For floating point data.
      double loc_dbl[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      double min_dbl[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      double max_dbl[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      double sum_dbl[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      std::string names[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      int k = 0;
      names[k] = "number_of_mapped_boxes (N)";
      loc_dbl[k++] = number_of_mapped_boxes;
      names[k] = "has_mapped_box";
      loc_dbl[k++] = has_mapped_box;
      names[k] = "smallest_aspect";
      loc_dbl[k++] = smallest_aspect;
      names[k] = "largest_aspect";
      loc_dbl[k++] = largest_aspect;
      names[k] = "sum_aspect";
      loc_dbl[k++] = sum_aspect;
      names[k] = "sum_normsurfarea";
      loc_dbl[k++] = sum_normsurfarea;
      names[k] = "sum_surfarea";
      loc_dbl[k++] = sum_surfarea;
      names[k] = "smallest_dim";
      loc_dbl[k++] = smallest_dim;
      names[k] = "largest_dim";
      loc_dbl[k++] = largest_dim;
      TBOX_ASSERT(k < MAPPED_BOX_LEVEL_NUMBER_OF_STATS);
      for (int i = 0; i < k; ++i) max_dbl[i] = sum_dbl[i] = min_dbl[i] =
                  loc_dbl[i];
      int rank_of_min[MAPPED_BOX_LEVEL_NUMBER_OF_STATS],
          rank_of_max[MAPPED_BOX_LEVEL_NUMBER_OF_STATS];
      if (d_mpi.getSize() > 1) {
         d_mpi.AllReduce(min_dbl, k, MPI_MINLOC, rank_of_min);
         d_mpi.AllReduce(max_dbl, k, MPI_MAXLOC, rank_of_max);
         d_mpi.AllReduce(sum_dbl, k, MPI_SUM);
      } else {
         for (int i = 0; i < k; ++i) {
            rank_of_min[i] = rank_of_max[i] = 0;
         }
      }

      co.unsetf(std::ios::fixed | std::ios::scientific);
      co.precision(3);

      co << border
      <<
      "                               local        min               max           sum/N    sum/P\n";
      for (int i = 0; i < k; ++i) {
         co << border << std::setw(27) << std::left << names[i]
         << ' ' << std::setw(8) << std::right << loc_dbl[i]
         << ' ' << std::setw(8) << std::right << min_dbl[i] << " @ "
         << std::setw(6) << std::left << rank_of_min[i]
         << ' ' << std::setw(8) << std::right << max_dbl[i] << " @ "
         << std::setw(6) << std::left << rank_of_max[i]
         << ' ' << std::setw(8) << std::right << sum_dbl[i] / getGlobalNumberOfBoxes()
         << ' ' << std::setw(8) << std::right << sum_dbl[i] / d_mpi.getSize()
         << '\n';
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
   t_acquire_remote_mapped_boxes = tbox::TimerManager::getManager()->
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
   t_acquire_remote_mapped_boxes.setNull();
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
