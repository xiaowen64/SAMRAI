/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of edges incident from a mapped_box_level of a distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_Connector_C
#define included_hier_Connector_C

#include "SAMRAI/hier/Connector.h"

#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/Connector.I"
#endif

#include <algorithm>
//#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

static const std::string dbgbord;

namespace SAMRAI {
namespace hier {

const int Connector::HIER_CONNECTOR_VERSION = 0;

tbox::Pointer<tbox::Timer> Connector::t_acquire_remote_relationships;
tbox::Pointer<tbox::Timer> Connector::t_cache_global_reduced_data;

tbox::StartupShutdownManager::Handler
Connector::s_initialize_finalize_handler(
   Connector::initializeCallback,
   0,
   0,
   Connector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector():
   d_base_handle(),
   d_head_handle(),
   d_base_width(tbox::Dimension::getInvalidDimension()),
   d_ratio(tbox::Dimension::getInvalidDimension()),
   d_ratio_is_exact(false),
   d_head_coarser(false),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_parallel_state(BoxLevel::DISTRIBUTED),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_connector_type(UNKNOWN)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const Connector& other):
   tbox::DescribedClass(),
   d_base_handle(other.d_base_handle),
   d_head_handle(other.d_head_handle),
   d_base_width(other.d_base_width),
   d_ratio(other.d_ratio),
   d_ratio_is_exact(other.d_ratio_is_exact),
   d_head_coarser(other.d_head_coarser),
   d_relationships(other.d_relationships),
   d_global_relationships(other.d_global_relationships),
   d_mpi(other.d_mpi),
   d_parallel_state(other.d_parallel_state),
   d_finalized(other.d_finalized),
   d_global_number_of_neighbor_sets(other.d_global_number_of_neighbor_sets),
   d_global_number_of_relationships(other.d_global_number_of_relationships),
   d_global_data_up_to_date(other.d_global_data_up_to_date),
   d_connector_type(other.d_connector_type)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const BoxLevel& base_mapped_box_level,
   const BoxLevel& head_mapped_box_level,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state):
   d_base_width(base_width.getDim(), 0),
   d_ratio(base_width.getDim(), 0),
   d_head_coarser(false),
   d_mpi(base_mapped_box_level.getMPI()),
   d_parallel_state(parallel_state),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(true),
   d_connector_type(UNKNOWN)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(base_mapped_box_level,
      head_mapped_box_level,
      base_width);

   setBase(base_mapped_box_level);
   setHead(head_mapped_box_level);
   setWidth(base_width, true);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::~Connector()
{
   clear();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

const NeighborhoodSet
& Connector::getGlobalNeighborhoodSets() const
{
   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Global connectivity unavailable in DISTRIBUTED state.");
   }
   return d_global_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::insertNeighbors(
   const NeighborSet& neighbors,
   const BoxId& mapped_box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != getBase().getMPI().getRank()) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(mapped_box_id)) {
      TBOX_ERROR(
         "Exiting due to above reported error."
         << "Connector::insertNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   const tbox::Dimension& dim = d_ratio.getDim();
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships[mapped_box_id].insert(neighbors.begin(),
         neighbors.end());
   }
   if (mapped_box_id.getOwnerRank() == getBase().getMPI().getRank()) {
      d_relationships[mapped_box_id].insert(neighbors.begin(),
         neighbors.end());
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::eraseNeighbor(
   const Box& neighbor,
   const BoxId& mapped_box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != getBase().getMPI().getRank()) {
      TBOX_ERROR("Connector::eraseNeighbor error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(mapped_box_id)) {
      TBOX_ERROR(
         "Connector::eraseNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      NeighborhoodSet::iterator mi(d_global_relationships.find(mapped_box_id));
      if (mi != d_global_relationships.end()) {
         mi->second.erase(neighbor);
      }
   }
   if (mapped_box_id.getOwnerRank() == getBase().getMPI().getRank()) {
      NeighborhoodSet::iterator mi(d_relationships.find(mapped_box_id));
      if (mi != d_relationships.end()) {
         mi->second.erase(neighbor);
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::shrinkWidth(const IntVector& new_width)
{
   if (!(new_width <= getConnectorWidth())) {
      TBOX_ERROR("Connector::shrinkWidth: new ghost cell\n"
         << "width " << new_width << " involves an\n"
         << "enlargement of the current cell width "
         << getConnectorWidth());
   }
   else if (new_width == getConnectorWidth()) {
      // This is a no-op.
      return;
   }

   // Have not yet written this for GLOBALIZED mode.
   TBOX_ASSERT(getParallelState() == BoxLevel::DISTRIBUTED);

   /*
    * Remove overlaps that disappeared given the new GCW.
    * Swap out the overlaps, modify them then swap them back in.
    */

   const bool head_coarser = getHeadCoarserFlag();
   const bool base_coarser = !getHeadCoarserFlag() &&
      getBase().getRefinementRatio() != getHead().getRefinementRatio();

   const tbox::ConstPointer<GridGeometry>& grid_geom(getBase().getGridGeometry());

   for (NeighborhoodSet::iterator ei = d_relationships.begin();
        ei != d_relationships.end(); ++ei) {
      const BoxId& mapped_box_id = ei->first;
      NeighborSet& nabrs = ei->second;
      const Box& mapped_box = *getBase().getBoxStrict(
            mapped_box_id);
      Box mapped_box_box = mapped_box;
      mapped_box_box.grow(new_width);
      if (base_coarser) {
         mapped_box_box.refine(getRatio());
      }
      for (BoxContainer::Iterator na = nabrs.begin();
           na != nabrs.end(); /* incremented in loop */) {
         const Box& nabr = *na;
         Box nabr_box = nabr;
         if (nabr.getBlockId() != mapped_box.getBlockId()) {
            grid_geom->transformBox(nabr_box,
               getHead().getRefinementRatio(),
               mapped_box.getBlockId(),
               nabr.getBlockId());
         }
         if (head_coarser) {
            nabr_box.refine(getRatio());
         }
         if (!mapped_box_box.intersects(nabr_box)) {
            nabrs.erase(na++);
         } else {
            ++na;
         }
      }
   }

   d_base_width = new_width;
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::removePeriodicRelationships()
{
   d_relationships.removePeriodicNeighbors();
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships.removePeriodicNeighbors();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::removePeriodicLocalNeighbors()
{
   for (NeighborhoodSet::iterator ei = d_relationships.begin();
        ei != d_relationships.end(); ++ei) {
      d_relationships[ei->first].removePeriodicImageBoxes();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool Connector::hasPeriodicLocalNeighborhoodRoots() const
{
   bool result = false;
   for (NeighborhoodSet::const_iterator ei = d_relationships.begin();
        ei != d_relationships.end(); ++ei) {
      if (ei->first.getPeriodicId().getPeriodicValue() != 0) {
         result = true;
         break;
      }
   }
   return result;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::acquireRemoteNeighborhoods()
{
   tbox::SAMRAI_MPI mpi(getMPI());
   if (mpi.getSize() == 1) {
      // In single-proc mode, we already have all the relationships already.
      d_global_relationships = d_relationships;
      return;
   }

   t_acquire_remote_relationships->start();

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack relationships from all mapped_box_level relationship sets into a single message.
    * Note that each mapped_box_level relationship set object packs the size of its
    * sub-message into send_mesg.
    */
   acquireRemoteNeighborhoods_pack(send_mesg,
      static_cast<int>(send_mesg.size()));
   int send_mesg_size = static_cast<int>(send_mesg.size());

   /*
    * Send and receive the data.
    */

   std::vector<int> recv_mesg_size(getMPI().getSize());
   mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(getMPI().getSize());
   int totl_size = 0;
   for (int n = 0; n < getMPI().getSize(); ++n) {
      proc_offset[n] = totl_size;
      totl_size += recv_mesg_size[n];
   }
   recv_mesg.resize(totl_size, BAD_INT);
   mpi.Allgatherv(&send_mesg[0],
      send_mesg_size,
      MPI_INT,
      &recv_mesg[0],
      &recv_mesg_size[0],
      &proc_offset[0],
      MPI_INT);

   /*
    * Extract relationship info received from other processors.
    */
   acquireRemoteNeighborhoods_unpack(recv_mesg, proc_offset);

   t_acquire_remote_relationships->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::acquireRemoteNeighborhoods_pack(
   std::vector<int>& send_mesg,
   int offset) const
{
   const tbox::Dimension dim(getBase().getDim());

   (void)offset;
   /*
    * relationship acquisition is done during globalization.
    * Thus, do not rely on current value of d_parallel_state.
    */

   /*
    * Pack relationship info from d_relationships into send_mesg,
    * starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of local mapped_boxes with neighbor lists (no info to send
    *     for those mapped_boxes without neighbor lists)
    *   - For each local mapped_box,
    *     - mapped_box's local index
    *     - number of neighbors
    *     - neighbors
    */
   const int num_mapped_boxes = static_cast<int>(getLocalNumberOfNeighborSets());
   int num_nabrs = 0;
   for (NeighborhoodSet::const_iterator ci = d_relationships.begin(); ci != d_relationships.end();
        ++ci) {
      num_nabrs += static_cast<int>((*ci).second.size());
   }
   const int mesg_size =
      1             /* number of local mapped_boxes with neighbor lists */
      + 3 * num_mapped_boxes /* local index, block id, and neighbor list size of each mapped_box */
      + num_nabrs * Box::commBufferSize(dim) /* neighbors */
   ;

   send_mesg.resize(mesg_size, BAD_INT);
   send_mesg[0] = num_mapped_boxes;
   int imesg = 1;

   for (NeighborhoodSet::const_iterator ci = d_relationships.begin(); ci != d_relationships.end();
        ++ci) {

      const BoxId& mapped_box_id = (*ci).first;
      const NeighborSet& nabrs = (*ci).second;

      send_mesg[imesg++] = mapped_box_id.getLocalId().getValue();
      send_mesg[imesg++] = mapped_box_id.getBlockId().getBlockValue();
      send_mesg[imesg++] = static_cast<int>(nabrs.size());

      for (NeighborSet::ConstIterator ni = nabrs.begin();
           ni != nabrs.end(); ++ni) {
         (*ni).putToIntBuffer(&send_mesg[imesg]);
         imesg += Box::commBufferSize(dim);
      }

   }

   TBOX_ASSERT(imesg == mesg_size);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::acquireRemoteNeighborhoods_unpack(
   const std::vector<int>& recv_mesg,
   const std::vector<int>& proc_offset)
{
   const tbox::Dimension dim(getBase().getDim());

   /*
    * Unpack relationship info from recv_mesg into d_relationships.
    */
   int mapped_box_com_buf_size = Box::commBufferSize(dim);

   d_global_relationships = d_relationships;

   for (int n = 0; n < getMPI().getSize(); ++n) {

      if (n != getMPI().getRank()) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int num_nabr_lists = *(ptr++);

         for (int i = 0; i < num_nabr_lists; ++i) {

            const LocalId local_id(*ptr++);
            const BlockId block_id(*ptr++);
            const int num_nabrs = (*ptr++);
            const BoxId mapped_box_id(local_id, n, block_id);

            NeighborSet& nabrs = d_global_relationships[mapped_box_id];
            Box nabr(dim);
            for (int nn = 0; nn < num_nabrs; ++nn) {
               nabr.getFromIntBuffer(ptr);
               nabrs.insert(nabrs.end(), nabr);
               ptr += mapped_box_com_buf_size;
            }

         }

      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::setParallelState(
   const BoxLevel::ParallelState parallel_state)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isFinalized()) {
      TBOX_ERROR(
         "Connector::setParallelState: Cannot change the parallel state of\n"
         << "an unfinalized Connector.  See Connector::finalizeContext()");
   }
#endif
   if (parallel_state != BoxLevel::DISTRIBUTED && parallel_state !=
       BoxLevel::GLOBALIZED) {
      TBOX_ERROR("Connector::setParallelState: Invalid distribution state: "
         << parallel_state << "\n");
   }

   if (d_parallel_state == BoxLevel::DISTRIBUTED && parallel_state ==
       BoxLevel::GLOBALIZED) {
      acquireRemoteNeighborhoods();
   } else if (d_parallel_state == BoxLevel::GLOBALIZED && parallel_state ==
              BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::finalizeContext()
{
   TBOX_ASSERT(!d_base_handle.isNull());
   TBOX_ASSERT(!d_head_handle.isNull());
   TBOX_ASSERT(d_base_width.getDim().isValid());

   const BoxLevel& base = d_base_handle->getBoxLevel();
   const BoxLevel& head = d_head_handle->getBoxLevel();
   const IntVector& baseRefinementRatio = base.getRefinementRatio();
   const IntVector& headRefinementRatio = head.getRefinementRatio();

   if (base.getGridGeometry() != head.getGridGeometry()) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Connector must be finalized with\n"
         << "BoxLevels using the same GridGeometry.");
   }
   if (!(baseRefinementRatio >= headRefinementRatio ||
         baseRefinementRatio <= headRefinementRatio)) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Refinement ratio between base and head box_levels\n"
         << "cannot be mixed (bigger in some dimension and\n"
         << "smaller in others).\n"
         << "Input base ratio = " << baseRefinementRatio
         << "\n"
         << "Input head ratio = " << headRefinementRatio
         << "\n");
   }
   if (d_parallel_state == BoxLevel::GLOBALIZED &&
       base.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR(
         "Connector::finalizeContext: base BoxLevel must be in\n"
         << "GLOBALIZED state for the Connector to be in\n"
         << "GLOBALIZED state.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   bool errf = false;
   for (NeighborhoodSet::const_iterator ci = d_relationships.begin();
        ci != d_relationships.end();
        ++ci) {
      if (!base.hasBox((*ci).first)) {
         const NeighborSet& nabrs = (*ci).second;
         tbox::perr << "\nConnector::finalizeContext: NeighborhoodSet "
                    << "provided for non-existent box " << ci->first
                    << "\n" << "Neighbors (" << nabrs.size() << "):\n";
         for (NeighborSet::ConstIterator na = nabrs.begin();
              na != nabrs.end(); ++na) {
            tbox::perr << (*na) << "\n";
         }
         errf = true;
      }
   }
   if (errf) {
      TBOX_ERROR(
         "Exiting due to errors."
         << "\nConnector::finalizeContext base box_level:\n"
         << base.format(dbgbord, 2));
   }
#endif
   computeRatioInfo(
      baseRefinementRatio,
      headRefinementRatio,
      d_ratio,
      d_head_coarser,
      d_ratio_is_exact);

   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   else {
      if (&d_relationships != &d_global_relationships) {
         d_global_relationships = d_relationships;
      }
   }

   // Erase remote relationships, if any, from d_relationships.
   if (!d_relationships.empty()) {
      if (d_relationships.begin()->first.getOwnerRank() != base.getMPI().getRank() ||
          d_relationships.rbegin()->first.getOwnerRank() != base.getMPI().getRank()) {
         NeighborhoodSet::Range range = d_relationships.findRanksRange(base.getMPI().getRank());
         if (range.first == range.second) {
            // No relationship belongs to local process.
            d_relationships.clear();
         }
         else {
            // Erase relationships belonging to remote process <
            // base.getMPI().getRank().
            d_relationships.erase(d_relationships.begin(), range.first);

            // Erase relationships belonging to remote process >
            // base.getMPI().getRank().
            d_relationships.erase(range.second, d_relationships.end());
         }
      }
   }
   d_finalized = true;
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::setBase(
   const BoxLevel& new_base,
   bool finalize_context)
{
   if (!new_base.isInitialized()) {
      TBOX_ERROR("Connector::setBase():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_base_handle = new_base.getBoxLevelHandle();
   d_mpi = new_base.getMPI();
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::setHead(
   const BoxLevel& new_head,
   bool finalize_context)
{
   if (!new_head.isInitialized()) {
      TBOX_ERROR("Connector::setHead():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_head_handle = new_head.getBoxLevelHandle();
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::setWidth(
   const IntVector& new_width,
   bool finalize_context)
{
   if (!(new_width >= IntVector::getZero(new_width.getDim()))) {
      TBOX_ERROR("Connector::setWidth():\n"
         << "Invalid ghost cell width: "
         << new_width << "\n");
   }
   d_finalized = false;
   d_base_width = new_width;
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::computeRatioInfo(
   const IntVector& baseRefinementRatio,
   const IntVector& headRefinementRatio,
   IntVector& ratio,
   bool& head_coarser,
   bool& ratio_is_exact)
{
   if (baseRefinementRatio <= headRefinementRatio) {
      ratio = headRefinementRatio / baseRefinementRatio;
      head_coarser = false;
      ratio_is_exact = (ratio * baseRefinementRatio) == headRefinementRatio;
   }
   else {
      ratio = baseRefinementRatio / headRefinementRatio;
      head_coarser = true;
      ratio_is_exact = (ratio * headRefinementRatio) == baseRefinementRatio;
   }
   if (baseRefinementRatio * headRefinementRatio <
       IntVector::getZero(baseRefinementRatio.getDim())) {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -headRefinementRatio * baseRefinementRatio;
      ratio_is_exact = true;
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::clear()
{
   if ( !d_base_handle.isNull() ) {
      d_relationships.clear();
      d_global_relationships.clear();
      d_mpi.setCommunicator(tbox::SAMRAI_MPI::commNull);
      d_base_handle.setNull();
      d_head_handle.setNull();
      d_base_width(0) = d_ratio(0) = 0;
      d_parallel_state = BoxLevel::DISTRIBUTED;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initializeToLocalTranspose(
   const Connector& connector)
{
   const IntVector my_gcw = convertHeadWidthToBase(
         connector.getHead().getRefinementRatio(),
         connector.getBase().getRefinementRatio(),
         connector.getConnectorWidth());

   clearNeighborhoods();
   setBase(connector.d_head_handle->getBoxLevel());
   setHead(connector.d_base_handle->getBoxLevel());
   setWidth(my_gcw, true);
   TBOX_ASSERT(isTransposeOf(connector));

   const tbox::Dimension dim(my_gcw.getDim());
   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   for (ConstNeighborhoodIterator ci = connector.begin();
        ci != connector.end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const BoxContainer::ConstIterator ni = getHead().getBox(mapped_box_id);
      if (ni == getHead().getBoxes().end()) {
         TBOX_ERROR(
            "Connector::initializeToLocalTranspose: mapped_box index\n"
            << mapped_box_id
            << " not found in local part of head mapped_box_level.\n"
            << "This means that the incoming Connector data was not a\n"
            << "self-consistent local mapping.\n");
      }
      const Box& my_head_mapped_box = *ni;

      for (ConstNeighborIterator na = connector.begin(ci);
           na != connector.end(ci); ++na) {
         const Box& my_base_mapped_box = *na;
         if (my_base_mapped_box.getOwnerRank() != getMPI().getRank()) {
            TBOX_ERROR(
               "Connector::initializeToLocalTranspose: base mapped_box "
               << my_head_mapped_box << "\n"
               << "has remote neighbor " << my_base_mapped_box
               << " which is disallowed.\n"
               << "Mapped_boxes must have only local neighbors in this method.");
         }
         if (my_base_mapped_box.isPeriodicImage()) {
            Box my_shifted_head_mapped_box(
               my_head_mapped_box,
               shift_catalog->getOppositeShiftNumber(
                  my_base_mapped_box.getPeriodicId()),
               getHead().getRefinementRatio());
            if (getHead().hasBox(my_shifted_head_mapped_box)) {
               BoxId base_non_per_id(
                  my_base_mapped_box.getGlobalId(),
                  my_base_mapped_box.getBlockId(),
                  PeriodicId::zero());
               d_relationships[base_non_per_id].insert(
                  my_shifted_head_mapped_box);
            }
         } else {
            d_relationships[my_base_mapped_box.getId()].insert(my_head_mapped_box);
         }
      }

   }

   if (0) {
      tbox::perr << "end of initializeToLocalTranspose:\n"
                 << "base:\n" << getBase().format("BASE->", 3)
                 << "head:\n" << getHead().format("HEAD->", 3)
                 << "this:\n" << format("THIS->", 3)
                 << "r:\n" << connector.format("RRRR->", 3)
                 << "Checking this transpose correctness:" << std::endl;
      assertTransposeCorrectness(connector, false);
      tbox::perr << "Checking r's transpose correctness:" << std::endl;
      connector.assertTransposeCorrectness(*this, false);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::eraseEmptyNeighborSets()
{
   for (NeighborhoodSet::iterator ei = d_relationships.begin();
        ei != d_relationships.end(); ) {
      if ((*ei).second.isEmpty()) {
         d_relationships.erase(ei++);
      } else {
         ++ei;
      }
   }
   d_global_data_up_to_date = false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
int
Connector::numLocalEmptyNeighborhoods() const
{
   int ct = 0;
   for (NeighborhoodSet::const_iterator itr = d_relationships.begin();
        itr != d_relationships.end(); ++itr) {
      if ((*itr).second.isEmpty()) {
         ++ct;
      }
   }
   return ct;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

bool
Connector::isTransposeOf(
   const Connector& other) const
{
   bool rval = false;
   if (d_base_handle == other.d_head_handle &&
       d_head_handle == other.d_base_handle) {
      if (d_head_coarser) {
         IntVector transpose_base_width = convertHeadWidthToBase(
               getHead().getRefinementRatio(),
               getBase().getRefinementRatio(),
               d_base_width);
         rval = other.d_base_width == transpose_base_width;
      } else {
         IntVector transpose_base_width = convertHeadWidthToBase(
               other.getHead().getRefinementRatio(),
               other.getBase().getRefinementRatio(),
               other.d_base_width);
         rval = d_base_width == transpose_base_width;
      }
   }
   return rval;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

bool Connector::isLocal() const
{
   for (NeighborhoodSet::const_iterator ei = d_relationships.begin(); ei != d_relationships.end();
        ++ei) {
      const BoxContainer& nabrs = ei->second;
      for (BoxContainer::ConstIterator na = nabrs.begin();
           na != nabrs.end();
           ++na) {
         if (na->getOwnerRank() != getMPI().getRank()) {
            return false;
         }
      }
   }
   return true;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t Connector::getLocalNumberOfNeighborSets() const
{
   return d_relationships.size();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t Connector::getLocalNumberOfRelationships() const
{
   size_t local_number_of_relationships = 0;
   for (NeighborhoodSet::const_iterator ei(d_relationships.begin());
        ei != d_relationships.end();
        ++ei) {
      local_number_of_relationships += static_cast<int>(ei->second.size());
   }
   return local_number_of_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
int Connector::getGlobalNumberOfNeighborSets() const
{
   TBOX_ASSERT(isFinalized());

   cacheGlobalReducedData();
   return d_global_number_of_neighbor_sets;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
int Connector::getGlobalNumberOfRelationships() const
{
   TBOX_ASSERT(isFinalized());

   cacheGlobalReducedData();
   return d_global_number_of_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isFinalized());

   if (d_global_data_up_to_date) {
      return;
   }

   t_cache_global_reduced_data->barrierAndStart();

   tbox::SAMRAI_MPI mpi(getMPI());

   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_number_of_relationships = 0;
      for (NeighborhoodSet::const_iterator ei(d_global_relationships.begin());
           ei != d_global_relationships.end();
           ++ei) {
         d_global_number_of_relationships +=
            static_cast<int>(ei->second.size());
      }
      d_global_number_of_neighbor_sets =
         static_cast<int>(d_global_relationships.size());
   } else {
      if (mpi.getSize() > 1) {
         int tmpa[2], tmpb[2];
         tmpa[0] = getLocalNumberOfNeighborSets();
         tmpa[1] = getLocalNumberOfRelationships();

         TBOX_ASSERT(tmpa[0] >= 0);
         TBOX_ASSERT(tmpa[0] >= 0);

         mpi.Allreduce(tmpa,
            tmpb,                        // Better to use MPI_IN_PLACE, but not some MPI's do not support.
            2,
            MPI_INT,
            MPI_SUM);
         d_global_number_of_neighbor_sets = tmpb[0];
         d_global_number_of_relationships = tmpb[1];
      } else {
         d_global_number_of_neighbor_sets = getLocalNumberOfNeighborSets();
         d_global_number_of_relationships = getLocalNumberOfRelationships();
      }

      TBOX_ASSERT(d_global_number_of_neighbor_sets >= 0);
      TBOX_ASSERT(d_global_number_of_relationships >= 0);
   }

   d_global_data_up_to_date = true;

   t_cache_global_reduced_data->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

IntVector Connector::convertHeadWidthToBase(
   const IntVector& base_refinement_ratio,
   const IntVector& head_refinement_ratio,
   const IntVector& head_gcw)
{
   if (!(base_refinement_ratio >= head_refinement_ratio ||
         base_refinement_ratio <= head_refinement_ratio)) {
      TBOX_ERROR("Connector::convertHeadWidthToBase:\n"
         << "head mapped_box_level must be either\n"
         << "finer or coarser than base.\n"
         << "Combined refinement and coarsening not allowed.");
   }

   tbox::Dimension dim(head_refinement_ratio.getDim());

   IntVector ratio(dim); // Ratio between head and base.

   if (head_refinement_ratio * base_refinement_ratio >
       IntVector::getZero(dim)) {
      // Same signs for both ratios -> simple to compute head-base ratio.
      if (base_refinement_ratio >= head_refinement_ratio) {
         ratio = base_refinement_ratio / head_refinement_ratio;
      } else {
         ratio = head_refinement_ratio / base_refinement_ratio;
      }
   } else {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -base_refinement_ratio * head_refinement_ratio;
   }
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));

   const IntVector base_width =
      (base_refinement_ratio >= head_refinement_ratio) ?
      (head_gcw * ratio) : IntVector::ceilingDivide(head_gcw, ratio);

   return base_width;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int detail_depth) const
{
   if (detail_depth < 0) {
      return;
   }

   if (!isFinalized()) {
      os << border << "Unfinalized.\n";
      return;
   }
   bool head_coarser = d_head_coarser;
   const IntVector head_gcw =
      convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         d_base_width);
   os << border << "Parallel state     : "
      << (getParallelState() == BoxLevel::DISTRIBUTED ? "DIST" : "GLOB")
      << '\n'
      << border << "Rank,nproc         : " << getMPI().getRank() << ", " << getMPI().getSize() << '\n'
      << border << "Base,head objects  :"
      << " ("
      << (d_base_handle == d_head_handle ? "same" : "different") << ") "
      << (void *)&d_base_handle->getBoxLevel() << ", "
      << (void *)&d_head_handle->getBoxLevel() << "\n"
      << border << "Base,head,/ ratios : "
      << getBase().getRefinementRatio() << ", "
      << getHead().getRefinementRatio() << ", "
      << d_ratio << (d_head_coarser ? " (head coarser)" : "") << '\n'
      << border << "Base,head widths   : " << d_base_width << ", "
      << head_gcw << '\n'
      << border << "Box count    : " << getBase().getLocalNumberOfBoxes()
      << " (" << getLocalNumberOfNeighborSets() << " with neighbor lists)\n"
   ;
   if (detail_depth > 0) {
      os << border << "Mapped_boxes with neighbors:\n";
      for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         BoxContainer::ConstIterator ni = getBase().getBox(ei->first);
         if (ni != getBase().getBoxes().end()) {
            os << border << "  "
               << (*ni) << "_"
               << (*ni).numberCells() << '\n';
         } else {
            os << border << "  #"
               << (*ei).first
               << ": INVALID DATA WARNING: nonexistent mapped_box index\n";
            TBOX_WARNING("Inconsistent data!!!\n"
               << "Neighbor data found for mapped_box "
               << (*ei).first << " but there is no such mapped_box!\n");
         }
         Box ghost_box = (*ni);
         ghost_box.grow(d_base_width);
         os << border << "    Neighbors (" << numLocalNeighbors(ei->first) << "):"
            << ((detail_depth > 1) ? "\n" : " ...\n");
         if (detail_depth > 1) {
            for (ConstNeighborIterator i_nabr = begin(ei);
                 i_nabr != end(ei); ++i_nabr) {
               Box ovlap = *i_nabr;
               if (ni->getBlockId() != i_nabr->getBlockId()) {
                  d_base_handle->getBoxLevel().getGridGeometry()->
                  transformBox(
                     ovlap,
                     d_head_handle->getBoxLevel().getRefinementRatio(),
                     ni->getBlockId(),
                     i_nabr->getBlockId());
               }
               if (head_coarser) ovlap.refine(d_ratio);
               else if (d_ratio != 1) ovlap.coarsen(d_ratio);
               ovlap = ovlap * ghost_box;
               os << border << "      "
                  << (*i_nabr) << "_"
                  << (*i_nabr).numberCells()
                  << "\tov" << ovlap << "_" << ovlap.numberCells() << '\n';
            }
         }
      }
   }
}

/*
 ***********************************************************************
 * Construct a Connector Outputter with formatting parameters.
 ***********************************************************************
 */

Connector::Outputter::Outputter(
   const Connector& connector,
   const std::string& border,
   int detail_depth):
   d_conn(connector),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a Connector according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
   std::ostream& os,
   const Connector::Outputter& format)
{
   format.d_conn.recursivePrint(os, format.d_border, format.d_detail_depth);
   return os;
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the Connector to a stream.
 ***********************************************************************
 */

Connector::Outputter Connector::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector *Connector::makeGlobalizedCopy(
   const Connector& other) const
{
   // Prevent wasteful accidental use when this method is not needed.
   TBOX_ASSERT(other.getParallelState() != BoxLevel::GLOBALIZED);

   Connector* copy = new Connector(other);
   copy->setParallelState(BoxLevel::GLOBALIZED);
   return copy;
}

/*
 ***********************************************************************
 * Run checkTransposeCorrectness and assert that no errors are found.
 ***********************************************************************
 */

void Connector::assertTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   size_t err_count =
      checkTransposeCorrectness(input_transpose, ignore_periodic_relationships);
   if (err_count) {
      TBOX_ERROR(
         "Connector::assertTransposeCorrectness:\n"
         << "Aborting with " << err_count << " transpose errors found:\n"
         << "this base:\n" << getBase().format("B:", 3)
         << "this head:\n" << getHead().format("H:", 3)
         << "this Connector:\n" << format("B->H:", 3)
         << "\ntranspose Connector:\n" << input_transpose.format("H->B:", 3));
   }
}

/*
 ***********************************************************************
 *
 * For every relationship in this, there should be reverse relationship in transpose.
 *
 * This method does not check whether the Connectors are defined to
 * form logical transposes (based on their widths and their base and
 * head mapped_box_levels).  For that, see isTransposeOf().
 *
 ***********************************************************************
 */

size_t Connector::checkTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   const tbox::Dimension dim(getBase().getDim());

   const Connector* transpose =
      (input_transpose.d_parallel_state == BoxLevel::GLOBALIZED) ?
      &input_transpose : makeGlobalizedCopy(input_transpose);

   const BoxLevel& head = getHead().getGlobalizedVersion();

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   /*
    * Check for extraneous relationships.
    * For every relationship in this, there should be reverse relationship in transpose.
    */
   Box shifted_mapped_box(dim);   // Shifted version of an unshifted Box.
   Box unshifted_mapped_box(dim); // Unhifted version of a shifted Box.

   size_t err_count = 0;


   const NeighborhoodSet& tran_relationships =
      transpose->getGlobalNeighborhoodSets();
   for (ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const Box& mapped_box = *getBase().getBox(mapped_box_id);

      size_t err_count_for_current_index = 0;

      for (ConstNeighborIterator ni = begin(ci); ni != end(ci); ++ni) {

         if (ignore_periodic_relationships && ni->isPeriodicImage()) {
            continue;
         }

         const Box& nabr = *ni;

         /*
          * Key for find in NeighborhoodSet must be non-periodic.
          */
         BoxId non_per_nabr_id(nabr.getGlobalId(),
                               nabr.getBlockId(),
                               PeriodicId::zero());

         const NeighborhoodSet::const_iterator cn =
            tran_relationships.find(non_per_nabr_id);

         if (cn == tran_relationships.end()) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local mapped_box " << mapped_box
            << " has relationship to " << nabr
            << " but " << nabr << " has no relationship container.\n";
            ++err_count_for_current_index;
            continue;
         }

         TBOX_ASSERT(cn->first == non_per_nabr_id);
         const NeighborSet& nabr_nabrs = cn->second;

         NeighborSet::ConstIterator nabr_ni(nabr_nabrs);

         if (nabr.isPeriodicImage()) {
            shifted_mapped_box.initialize(
               mapped_box,
               shift_catalog->getOppositeShiftNumber(nabr.getPeriodicId()),
               getBase().getRefinementRatio());
            nabr_ni = nabr_nabrs.find(shifted_mapped_box);
         } else {
            nabr_ni = nabr_nabrs.find(mapped_box);
         }

         if (nabr_ni == nabr_nabrs.end()) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local mapped_box " << mapped_box;
            if (nabr.isPeriodicImage()) {
               tbox::perr << " (shifted version " << shifted_mapped_box << ")";
            }
            tbox::perr << " has relationship to " << nabr << " but "
            << nabr << " does not have the reverse relationship.\n"
            ;
            tbox::perr << "Neighbors of " << nabr << " are:\n";
            for (NeighborSet::ConstIterator nj = nabr_nabrs.begin();
                 nj != nabr_nabrs.end(); ++nj) {
               tbox::perr << *nj << std::endl;
            }
            ++err_count_for_current_index;
            continue;
         }

      }

      if (err_count_for_current_index > 0) {
         tbox::perr << "Mapped_box " << mapped_box << " had "
         << err_count_for_current_index
         << " errors.  Neighbors are:\n";
         for (ConstNeighborIterator nj = begin(ci); nj != end(ci); ++nj) {
            tbox::perr << *nj << std::endl;
         }
         err_count += err_count_for_current_index;
      }

   }

   /*
    * Check for missing relationships:
    * Transpose should not contain any relationship that does not correspond to
    * one in this.
    */

   for (NeighborhoodSet::const_iterator ci = tran_relationships.begin();
        ci != tran_relationships.end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const NeighborSet& nabrs = ci->second;

      size_t err_count_for_current_index = 0;

      if (!head.hasBox(mapped_box_id)) {
         TBOX_ASSERT(head.hasBox(mapped_box_id));
      }
      const Box& head_mapped_box = *head.getBoxStrict(mapped_box_id);

      for (NeighborSet::ConstIterator na = nabrs.begin();
           na != nabrs.end(); ++na) {

         const Box nabr = *na;

         if (nabr.getOwnerRank() == getMPI().getRank()) {

            if (ignore_periodic_relationships && nabr.isPeriodicImage()) {
               continue;
            }

            if (!getBase().hasBox(nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box
               << " has neighbor " << nabr << "\n"
               << " but the neighbor does not exist "
               << "in the base mapped_box_level.\n";
               tbox::perr << "Neighbors of head mapped_box "
               << mapped_box_id << " are:\n";
               for (NeighborSet::ConstIterator nj = nabrs.begin();
                    nj != nabrs.end(); ++nj) {
                  tbox::perr << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box& base_mapped_box = *getBase().getBoxStrict(
                  nabr);

            /*
             * Non-periodic BoxId needed for NeighborhoodSet::find()
             */
            BoxId base_non_per_id(base_mapped_box.getGlobalId(),
                                  base_mapped_box.getBlockId(),
                                  PeriodicId::zero());
            NeighborhoodSet::const_iterator nabr_nabrs_ =
               d_relationships.find(base_non_per_id);

            if (nabr_nabrs_ == d_relationships.end()) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box << "\n"
               << " has base mapped_box "
               << base_mapped_box << " as a neighbor.\n"
               << "But " << base_mapped_box
               << " has no neighbor container.\n";
               tbox::perr << "Neighbors of head mapped_box " << BoxId(
                  mapped_box_id)
               << ":" << std::endl;
               for (NeighborSet::ConstIterator nj = nabrs.begin();
                    nj != nabrs.end(); ++nj) {
                  tbox::perr << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const NeighborSet& nabr_nabrs = nabr_nabrs_->second;

            const Box nabr_nabr(dim, mapped_box_id.getGlobalId(),
                                head_mapped_box.getBlockId(),
                                shift_catalog->getOppositeShiftNumber(
                                   base_mapped_box.getPeriodicId()));

            NeighborSet::ConstIterator found_nabr_ =
               nabr_nabrs.find(nabr_nabr);

            if (found_nabr_ == nabr_nabrs.end()) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box << "\n"
               << " has base mapped_box " << base_mapped_box
               << " as a neighbor.\n"
               << "But base mapped_box " << base_mapped_box
               << " does not have a mapped_box indexed "
               << nabr_nabr.getId()
               << " in its neighbor list." << std::endl;
               tbox::perr << "Neighbors of head mapped_box " << nabr_nabr.getId()
               << ":" << std::endl;
               for (NeighborSet::ConstIterator
                    nj = nabrs.begin(); nj != nabrs.end(); ++nj) {
                  tbox::perr << *nj << std::endl;
               }
               tbox::perr << "Neighbors of base mapped_box ";
               if (nabr.isPeriodicImage()) {
                  unshifted_mapped_box.initialize(
                     nabr,
                     shift_catalog->getZeroShiftNumber(),
                     getBase().getRefinementRatio());
                  tbox::perr << unshifted_mapped_box;
               }
               tbox::perr << ":" << std::endl;
               for (NeighborSet::ConstIterator nj = nabr_nabrs.begin();
                    nj != nabr_nabrs.end(); ++nj) {
                  tbox::perr << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

         }

      }

      if (err_count_for_current_index > 0) {
         err_count += err_count_for_current_index;
      }

   }

   if (transpose != &input_transpose) {
      delete transpose;
   }

   return err_count;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t Connector::checkConsistencyWithBase() const
{
   size_t num_errors = 0;
   for (ConstNeighborhoodIterator i_relationships = begin();
        i_relationships != end(); ++i_relationships) {
      const BoxId& mapped_box_id = (*i_relationships).first;
      if (!getBase().hasBox(mapped_box_id)) {
         ++num_errors;
         tbox::plog << "ERROR->"
         << "Connector::assertConsistencyWithBase: Neighbor data given "
         << "\nfor mapped_box " << mapped_box_id
         << " but the mapped_box does not exist.\n";
      }
   }
   return num_errors;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::assertConsistencyWithBase() const
{
   if (checkConsistencyWithBase() > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithBase() found inconsistencies.\n"
         << "Base mapped box level:\n" << getBase().format("ERROR->", 2));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::computeNeighborhoodDifferences(
   Connector& left_minus_right,
   const Connector& left,
   const Connector& right)
{
//TODO:  Commenting out until there's a way to do set_difference for
// BoxContainer
#if 1
   if (0) {
      tbox::plog << "Computing relationship differences, a:\n" << left.format(dbgbord, 3)
      << "Computing relationship differences, b:\n" << right.format(dbgbord, 3);
   }
   left_minus_right.clearNeighborhoods();
   left_minus_right.d_parallel_state = left.getParallelState();
   left_minus_right.setBase(left.d_base_handle->getBoxLevel());
   left_minus_right.setHead(left.d_head_handle->getBoxLevel());
   left_minus_right.setWidth(left.d_base_width, true);
   NeighborhoodSet& drelationships = left_minus_right.d_relationships;

   for (ConstNeighborhoodIterator ai = left.begin(); ai != left.end(); ++ai) {

      const BoxId& mapped_box_id = ai->first;
      const NeighborSet& anabrs = ai->second;

      ConstNeighborhoodIterator bi = right.findLocal(mapped_box_id);
      if (bi != right.end()) {
         const NeighborSet& bnabrs = bi->second;
         // Remove bi from ai.  Put results in a_minus_b.
         NeighborSet& diff = drelationships[mapped_box_id];

         //equivalent of stl set_difference
         BoxContainer::ConstIterator na = anabrs.begin();
         BoxContainer::ConstIterator nb = bnabrs.begin();
         while (na != anabrs.end() && nb != bnabrs.end()) {
            if ((*na).getId() < (*nb).getId()) {
               diff.insert(diff.end(), (*na));
               ++na;
            } else if ((*nb).getId() < (*na).getId()) {
               ++nb;
            } else {
               ++na;
               ++nb;
            } 
         }
             

#if 0
         std::insert_iterator<NeighborSet> ii(diff, diff.begin());
         set_difference(anabrs.begin(),
            anabrs.end(),
            bnabrs.begin(),
            bnabrs.end(),
            ii, Box::id_less());
#endif
         if (diff.isEmpty()) {
            drelationships.erase(mapped_box_id);
         }
      } else if (!anabrs.isEmpty()) {
         drelationships[mapped_box_id] = anabrs;
      }

   }
#endif
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::assertConsistencyWithHead() const
{
   const int number_of_inconsistencies = checkConsistencyWithHead();
   if (number_of_inconsistencies > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithHead() found inconsistencies.\n"
         << getBase().format("base-> ", 3)
         << getHead().format("head-> ", 3)
         << format("E-> ", 3));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t Connector::checkConsistencyWithHead() const
{
   const BoxLevel& head_mapped_box_level = getHead().getGlobalizedVersion();

   TBOX_ASSERT(head_mapped_box_level.getParallelState() ==
      BoxLevel::GLOBALIZED);

   const BoxContainer& head_mapped_boxes =
      head_mapped_box_level.getGlobalBoxes();

   size_t number_of_inconsistencies = 0;

   /*
    * For each neighbor in each neighbor list,
    * check that the neighbor is in the head_mapped_box_level.
    */

   for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {

      const BoxId& mapped_box_id = ei->first;

      for (ConstNeighborIterator na = begin(ei); na != end(ei); ++na) {

         const Box& nabr = *na;
         const Box unshifted_nabr(
            nabr, PeriodicId::zero(), head_mapped_box_level.getRefinementRatio());

         BoxContainer::ConstIterator na_in_head =
            head_mapped_boxes.find(unshifted_nabr);

         if (na_in_head == head_mapped_boxes.end()) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Neighbor list for mapped_box " << mapped_box_id << "\n"
            << "referenced nonexistent neighbor "
            << nabr << "\n";
            tbox::perr << "Neighbors of mapped_box " << mapped_box_id << ":\n";
            for (ConstNeighborIterator nb = begin(ei); nb != end(ei); ++nb) {
               tbox::perr << "    " << *nb << '\n';
            }
            ++number_of_inconsistencies;
            continue;
         }

         const Box& nabr_in_head = *na_in_head;
         if (!unshifted_nabr.isIdEqual(nabr_in_head) ||
             !unshifted_nabr.isSpatiallyEqual(nabr_in_head)) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Inconsistent mapped_box data at mapped_box "
            << mapped_box_id << "\n"
            << "Neighbor " << nabr << "(unshifted to "
            << unshifted_nabr << ") does not match "
            << "head mapped_box " << nabr_in_head
            << "\n";
            ++number_of_inconsistencies;
         }

      }
   }

   return number_of_inconsistencies;
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void Connector::getNeighborBoxes(
   const BoxId& mapped_box_id,
   BoxContainer& nbr_boxes) const
{
   const NeighborSet& nbr_mapped_boxes = getNeighborSet(mapped_box_id);
   for (NeighborSet::ConstIterator ni = nbr_mapped_boxes.begin();
        ni != nbr_mapped_boxes.end(); ni++) {
      nbr_boxes.pushBack(*ni);
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void Connector::setConnectorType(
   ConnectorType connector_type)
{
   d_connector_type = connector_type;
}

/*
 ***************************************************************************
 ***************************************************************************
 */
Connector::ConnectorType Connector::getConnectorType() const
{
   return d_connector_type;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initializeCallback()
{
   t_acquire_remote_relationships = tbox::TimerManager::getManager()->
      getTimer("hier::Connector::acquireRemoteNeighborhoods()");
   t_cache_global_reduced_data = tbox::TimerManager::getManager()->
      getTimer("hier::Connector::cacheGlobalReducedData()");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */

void Connector::finalizeCallback()
{
   t_acquire_remote_relationships.setNull();
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
