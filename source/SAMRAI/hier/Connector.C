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

tbox::Pointer<tbox::Timer> Connector::t_initialize_private;
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
   d_relationships(tbox::Dimension::getInvalidDimension()),
   d_global_relationships(tbox::Dimension::getInvalidDimension()),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_parallel_state(BoxLevel::DISTRIBUTED),
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
   d_relationships(base_width.getDim()),
   d_global_relationships(base_width.getDim()),
   d_mpi(base_mapped_box_level.getMPI()),
   d_parallel_state(BoxLevel::DISTRIBUTED),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(true),
   d_connector_type(UNKNOWN)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(base_mapped_box_level,
      head_mapped_box_level,
      base_width);

   initialize(base_mapped_box_level,
      head_mapped_box_level,
      base_width,
      parallel_state);
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
   TBOX_ASSERT(isInitialized());
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
      NeighborSet& global_nabrs =
         d_global_relationships.getNeighborSet(
            mapped_box_id, dim);
      global_nabrs.insert(neighbors.orderedBegin(),
                          neighbors.orderedEnd());
   }
   if (mapped_box_id.getOwnerRank() == getBase().getMPI().getRank()) {
      NeighborSet& nabrs = d_relationships.getNeighborSet(mapped_box_id, dim);
      nabrs.insert(neighbors.orderedBegin(),
                   neighbors.orderedEnd());
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
      for (BoxSet::OrderedConstIterator na = nabrs.orderedBegin();
           na != nabrs.orderedEnd(); /* incremented in loop */) {
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
   TBOX_ASSERT(isInitialized());
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
   TBOX_ASSERT(isInitialized());
   for (NeighborhoodSet::iterator ei = d_relationships.begin();
        ei != d_relationships.end(); ++ei) {
      d_relationships.getNeighborSet(ei->first, d_ratio.getDim()).removePeriodicImageBoxes();
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

   std::vector<int> recv_mesg_size(getBase().getMPI().getSize());
   mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(getBase().getMPI().getSize());
   int totl_size = 0;
   for (int n = 0; n < getBase().getMPI().getSize(); ++n) {
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

      for (NeighborSet::OrderedConstIterator ni = nabrs.orderedBegin();
           ni != nabrs.orderedEnd(); ++ni) {
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

   for (int n = 0; n < getBase().getMPI().getSize(); ++n) {

      if (n != getBase().getMPI().getRank()) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int num_nabr_lists = *(ptr++);

         for (int i = 0; i < num_nabr_lists; ++i) {

            const LocalId local_id(*ptr++);
            const BlockId block_id(*ptr++);
            const int num_nabrs = (*ptr++);
            const BoxId mapped_box_id(local_id, n, block_id);

            NeighborSet& nabrs =
               d_global_relationships.getNeighborSet(mapped_box_id, dim);
            Box nabr(dim);
            for (int nn = 0; nn < num_nabrs; ++nn) {
               nabr.getFromIntBuffer(ptr);
               nabrs.insert(nabrs.orderedEnd(), nabr);
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
   if (!isInitialized()) {
      TBOX_ERROR(
         "Connector::setParallelState: Cannot change the parallel state of\n"
         << "an uninitialized Connector.  See Connector::initialize()");
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

void Connector::initialize(
   const BoxLevel& base,
   const BoxLevel& head,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state,
   bool clear_relationships)
{
   if (clear_relationships) {
      d_relationships.clear();
   }
   initializePrivate(base, head, base_width, base.getRefinementRatio(),
      head.getRefinementRatio(), parallel_state);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initializePrivate(
   const BoxLevel& base,
   const BoxLevel& head,
   const IntVector& base_width,
   const IntVector& baseRefinementRatio,
   const IntVector& headRefinementRatio,
   const BoxLevel::ParallelState parallel_state)
{
   t_initialize_private->start();
   /*
    * Check inputs.
    */
   if (!base.isInitialized() ||
       !head.isInitialized()) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Connector may not be initialized with\n"
         << "an uninitialized BoxLevel.");
   }
   if (base.getGridGeometry() != head.getGridGeometry()) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Connector must be initialized with\n"
         << "BoxLevels using the same GridGeometry.");
   }
   if (!(base.getRefinementRatio() >=
         head.getRefinementRatio() ||
         base.getRefinementRatio() <=
         head.getRefinementRatio())) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Refinement ratio between base and head mapped_box_levels\n"
         << "cannot be mixed (bigger in some dimension and\n"
         << "smaller in others).\n"
         << "Input base ratio = " << base.getRefinementRatio()
         << "\n"
         << "Input head ratio = " << head.getRefinementRatio()
         << "\n");
   }
   if (!(base_width >= IntVector::getZero(base_width.getDim()))) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Invalid ghost cell width: "
         << base_width << "\n");
   }
   if (parallel_state == BoxLevel::GLOBALIZED &&
       base.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR(
         "Connector::initializePrivate: base BoxLevel must be in\n"
         << "GLOBALIZED state before initializing the Connector to\n"
         << "GLOBALIZED state.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   bool errf = false;
   for (NeighborhoodSet::const_iterator ci = d_relationships.begin();
        ci != d_relationships.end();
        ++ci) {
      if (!base.hasBox((*ci).first)) {
         const NeighborSet& nabrs = (*ci).second;
         tbox::perr << "\nConnector::initializePrivate: NeighborhoodSet "
                    << "provided for non-existent mapped_box " << ci->first
                    << "\n" << "Neighbors (" << nabrs.size() << "):\n";
         for (NeighborSet::OrderedConstIterator na = nabrs.orderedBegin();
              na != nabrs.orderedEnd(); ++na) {
            tbox::perr << (*na) << "\n";
         }
         errf = true;
      }
   }
   if (errf) {
      TBOX_ERROR(
         "Exiting due to errors."
         << "\nConnector::initializePrivate base mapped_box_level:\n"
         << base.format(dbgbord, 2)
         << "\nConnector::initializePrivate head mapped_box_level:\n"
         << head.format(dbgbord, 2));
   }
#endif

   d_base_handle = base.getBoxLevelHandle();
   d_head_handle = head.getBoxLevelHandle();
   d_mpi = base.getMPI();

   d_base_width = base_width;

   computeRatioInfo(baseRefinementRatio, headRefinementRatio,
      d_ratio, d_head_coarser, d_ratio_is_exact);

   if (parallel_state == BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   } else {
      if (&d_relationships != &d_global_relationships) {
         d_global_relationships = d_relationships;
      }
   }

   // Erase remote relationships, if any, from d_relationships.
   if (!d_relationships.empty()) {
      if (d_relationships.begin()->first.getOwnerRank() != getBase().getMPI().getRank() ||
          d_relationships.rbegin()->first.getOwnerRank() != getBase().getMPI().getRank()) {
         NeighborhoodSet::Range range = d_relationships.findRanksRange(getBase().getMPI().getRank());
         if (range.first == range.second) {
            // No relationship belongs to local process.
            d_relationships.clear();
         } else {
            // Erase relationships belonging to remote process < getBase().getMPI().getRank().
            d_relationships.erase(d_relationships.begin(), range.first);
            // Erase relationships belonging to remote process > getBase().getMPI().getRank().
            d_relationships.erase(range.second, d_relationships.end());
         }
      }
   }

   d_parallel_state = parallel_state;

   t_initialize_private->stop();
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

   initialize(
      connector.d_head_handle->getBoxLevel(),
      connector.d_base_handle->getBoxLevel(),
      my_gcw,
      BoxLevel::DISTRIBUTED);
   TBOX_ASSERT(isTransposeOf(connector));

   const tbox::Dimension dim(my_gcw.getDim());
   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   for (ConstNeighborhoodIterator ci = connector.begin();
        ci != connector.end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const BoxSet::OrderedConstIterator ni = getHead().getBox(mapped_box_id);
      if (ni == getHead().getBoxes().orderedEnd()) {
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
         if (my_base_mapped_box.getOwnerRank() != getBase().getMPI().getRank()) {
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
               d_relationships.insertNeighbor(base_non_per_id,
                                              my_shifted_head_mapped_box);
            }
         } else {
            d_relationships.insertNeighbor(my_base_mapped_box.getId(),
                                           my_head_mapped_box);
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
      const BoxSet& nabrs = ei->second;
      for (BoxSet::ConstIterator na = nabrs.begin();
           na != nabrs.end();
           ++na) {
         if ((*na).getOwnerRank() != getBase().getMPI().getRank()) {
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
   TBOX_ASSERT(isInitialized());
   return d_relationships.size();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t Connector::getLocalNumberOfRelationships() const
{
   TBOX_ASSERT(isInitialized());
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
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_neighbor_sets;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
int Connector::getGlobalNumberOfRelationships() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void Connector::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isInitialized());

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
   if (detail_depth < 0) return;

   if (!isInitialized()) {
      os << border << "Uninitialized.\n";
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
      << border << "Rank,nproc         : " << getBase().getMPI().getRank() << ", " << getBase().getMPI().getSize() << '\n'
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
         BoxSet::OrderedConstIterator ni = getBase().getBox(ei->first);
         if (ni != getBase().getBoxes().orderedEnd()) {
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
 * Write out some statistics on the relationships, including statistics
 * for judging mesh quality.
 ***********************************************************************
 */

void Connector::printNeighborStats(
   std::ostream& os,
   const std::string& border) const
{
   if (!isInitialized()) {
      os << "Uninitialized.\n";
      return;
   }

   tbox::SAMRAI_MPI mpi(getMPI());

   getBase().cacheGlobalReducedData();

   /*
    * Compute additional statistics.
    */
   std::set<int> owners; // Neighbor owners for the local processor.

   int sum_nabr_sets = 0;
   int sum_relationships = 0;
   int sum_local_relationships = 0;
   int sum_nonlocal_relationships = 0;
   int sum_cells = 0;
   int sum_ovlap_size = 0;
   int sum_local_ovlap_size = 0;
   int sum_nonlocal_ovlap_size = 0;
   int sum_owners = 0;
   int sum_local_owners = 0;
   int sum_nonlocal_owners = 0;

   int max_nabrs = 0;
   int min_nabrs = 9999;

   for (ConstNeighborhoodIterator ni = begin(); ni != end(); ++ni) {

      ++sum_nabr_sets;
      const Box& mapped_box = *getBase().getBoxStrict(ni->first);
      sum_cells += mapped_box.size();
      Box base_ghost_box = mapped_box;
      base_ghost_box.grow(d_base_width);

      int num_local_nbrs = numLocalNeighbors(ni->first);
      max_nabrs = tbox::MathUtilities<int>::Max(max_nabrs, num_local_nbrs);
      min_nabrs = tbox::MathUtilities<int>::Min(min_nabrs, num_local_nbrs);

      for (ConstNeighborIterator na = begin(ni); na != end(ni); ++na) {

         const Box& nabr = *na;

         Box head_box = nabr;
         if (d_head_coarser) {
            head_box.refine(d_ratio);
         } else {
            head_box.coarsen(d_ratio);
         }
         Box overlap_box = base_ghost_box * head_box;
         int overlap_size = overlap_box.size();

         ++sum_relationships;
         sum_ovlap_size += overlap_size;

         if (nabr.getOwnerRank() == getBase().getMPI().getRank()) {
            ++sum_local_relationships;
            sum_local_ovlap_size += overlap_size;
            sum_local_owners = 1;
         } else {
            ++sum_nonlocal_relationships;
            sum_nonlocal_ovlap_size += overlap_size;
            owners.insert(nabr.getOwnerRank());
         }

      }

   }

   sum_nonlocal_owners = static_cast<int>(owners.size());
   sum_owners = sum_local_owners + sum_nonlocal_owners;

   /*
    * Collect local statistics into arrays for global reduction.
    */
   const int number_of_stats(20);
   int loc_num[number_of_stats];
   int min_num[number_of_stats];
   int max_num[number_of_stats];
   int sum_num[number_of_stats];
   std::string names[number_of_stats];
   int k = 0;
   names[k] = "sum_nabr_sets";
   loc_num[k++] = sum_nabr_sets;
   names[k] = "sum_relationships";
   loc_num[k++] = sum_relationships;
   const int k_sum_relationships = k - 1;
   names[k] = "sum_local_relationships";
   loc_num[k++] = sum_local_relationships;
   const int k_sum_local_relationships = k - 1;
   names[k] = "sum_nonlocal_relationships";
   loc_num[k++] = sum_nonlocal_relationships;
   names[k] = "sum_cells";
   loc_num[k++] = sum_cells;
   names[k] = "sum_ovlap_size";
   loc_num[k++] = sum_ovlap_size;
   const int k_sum_ovlap_size = k - 1;
   names[k] = "sum_local_ovlap_size";
   loc_num[k++] = sum_local_ovlap_size;
   const int k_sum_local_ovlap_size = k - 1;
   names[k] = "sum_nonlocal_ovlap_size";
   loc_num[k++] = sum_nonlocal_ovlap_size;
   names[k] = "sum_owners";
   loc_num[k++] = sum_owners;
   names[k] = "sum_local_owners";
   loc_num[k++] = sum_local_owners;
   names[k] = "sum_nonlocal_owners";
   loc_num[k++] = sum_nonlocal_owners;
   names[k] = "min_nabrs";
   loc_num[k++] = min_nabrs;
   names[k] = "max_nabrs";
   loc_num[k++] = max_nabrs;
   TBOX_ASSERT(k < number_of_stats);
   for (int i = 0; i < k; ++i) max_num[i] = sum_num[i] = min_num[i] =
               loc_num[i];
   int rank_of_min[number_of_stats],
       rank_of_max[number_of_stats];
   if (getBase().getMPI().getSize() > 1) {
      mpi.AllReduce(min_num, k, MPI_MINLOC, rank_of_min);
      mpi.AllReduce(max_num, k, MPI_MAXLOC, rank_of_max);
      mpi.AllReduce(sum_num, k, MPI_SUM);
   } else {
      for (int i = 0; i < k; ++i) {
         rank_of_min[i] = rank_of_max[i] = 0;
      }
   }

   os.unsetf(std::ios::fixed | std::ios::scientific);
   os.precision(3);

   os << border
      <<
   "                               local        min               max           sum/N    sum/P\n";
   for (int i = 0; i < k; ++i) {
      os << border << std::setw(27) << std::left << names[i]
         << ' ' << std::setw(8) << std::right << loc_num[i]
         << ' ' << std::setw(8) << std::right << min_num[i] << " @ " << std::setw(
         6) << std::left << rank_of_min[i]
         << ' ' << std::setw(8) << std::right << max_num[i] << " @ " << std::setw(
         6) << std::left << rank_of_max[i]
         << ' ' << std::setw(8) << std::right
         << (getBase().getGlobalNumberOfBoxes() > 0 ? static_cast<float>(sum_num[i]
                                                                      / getBase().
                                                                      getGlobalNumberOfBoxes()) :
          0.0)
      << ' ' << std::setw(8) << std::right << static_cast<float>(sum_num[i] / getBase().getMPI().getSize())
      << '\n';
   }

   // Write statistics that cannot be expressed as average over P processors or N Boxes.
   os << border
      << "                               local        min               max             avg\n";
   double local_relationship_fraction = loc_num[k_sum_relationships] ==
      0 ? 0.0 : (double)loc_num[k_sum_local_relationships] / loc_num[k_sum_relationships];
   double min_local_relationship_fraction = local_relationship_fraction;
   double max_local_relationship_fraction = local_relationship_fraction;
   double sum_local_relationship_fraction = local_relationship_fraction;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min_local_relationship_fraction, 1, MPI_MINLOC, rank_of_min);
      mpi.AllReduce(&max_local_relationship_fraction, 1, MPI_MAXLOC, rank_of_max);
      mpi.AllReduce(&sum_local_relationship_fraction, 1, MPI_SUM);
   }
   os << border << std::setw(27) << std::left << "local_relationship_fraction"
      << ' ' << std::setw(8) << std::right << local_relationship_fraction
      << ' ' << std::setw(8) << std::right << min_local_relationship_fraction << " @ "
      << std::setw(6) << std::left << rank_of_min[0]
      << ' ' << std::setw(8) << std::right << max_local_relationship_fraction << " @ "
      << std::setw(6) << std::left << rank_of_max[0]
      << ' ' << std::setw(8) << std::right
      << (sum_num[k_sum_relationships] ==
       0 ? 0.0 : (double)sum_num[k_sum_local_relationships] / sum_num[k_sum_relationships])
      << '\n';
   double local_ovlap_fraction = loc_num[k_sum_ovlap_size] ==
      0 ? 0.0 : (double)loc_num[k_sum_local_ovlap_size]
      / loc_num[k_sum_ovlap_size];
   double min_local_ovlap_fraction = local_ovlap_fraction;
   double max_local_ovlap_fraction = local_ovlap_fraction;
   double sum_local_ovlap_fraction = local_ovlap_fraction;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min_local_ovlap_fraction, 1, MPI_MINLOC, rank_of_min);
      mpi.AllReduce(&max_local_ovlap_fraction, 1, MPI_MAXLOC, rank_of_max);
      mpi.AllReduce(&sum_local_ovlap_fraction, 1, MPI_SUM);
   }
   os << border << std::setw(27) << std::left << "local_ovlap_fraction"
      << ' ' << std::setw(8) << std::right << local_ovlap_fraction
      << ' ' << std::setw(8) << std::right << min_local_ovlap_fraction << " @ "
      << std::setw(6) << std::left << rank_of_min[0]
      << ' ' << std::setw(8) << std::right << max_local_ovlap_fraction << " @ "
      << std::setw(6) << std::left << rank_of_max[0]
      << ' ' << std::setw(8) << std::right
      << (sum_num[k_sum_ovlap_size] ==
       0 ? 0.0 : (double)sum_num[k_sum_local_ovlap_size]
       / sum_num[k_sum_ovlap_size])
   << '\n';
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

   for (ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const Box& mapped_box = *getBase().getBox(mapped_box_id);

      size_t err_count_for_current_index = 0;

      for (ConstNeighborIterator ni = begin(ci); ni != end(ci); ++ni) {

         if (ignore_periodic_relationships && ni->isPeriodicImage()) {
            continue;
         }

         const Box& nabr = *ni;

         const NeighborhoodSet& tran_relationships = transpose->getGlobalNeighborhoodSets();

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

         NeighborSet::OrderedConstIterator nabr_ni(nabr_nabrs);

         if (nabr.isPeriodicImage()) {
            shifted_mapped_box.initialize(
               mapped_box,
               shift_catalog->getOppositeShiftNumber(nabr.getPeriodicId()),
               getBase().getRefinementRatio());
            nabr_ni = nabr_nabrs.find(shifted_mapped_box);
         } else {
            nabr_ni = nabr_nabrs.find(mapped_box);
         }

         if (nabr_ni == nabr_nabrs.orderedEnd()) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local mapped_box " << mapped_box;
            if (nabr.isPeriodicImage()) {
               tbox::perr << " (shifted version " << shifted_mapped_box << ")";
            }
            tbox::perr << " has relationship to " << nabr << " but "
            << nabr << " does not have the reverse relationship.\n"
            ;
            tbox::perr << "Neighbors of " << nabr << " are:\n";
            for (NeighborSet::OrderedConstIterator nj = nabr_nabrs.orderedBegin();
                 nj != nabr_nabrs.orderedEnd(); ++nj) {
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

   const NeighborhoodSet& tran_relationships = transpose->getGlobalNeighborhoodSets();

   for (NeighborhoodSet::const_iterator ci = tran_relationships.begin();
        ci != tran_relationships.end(); ++ci) {

      const BoxId& mapped_box_id = ci->first;
      const NeighborSet& nabrs = ci->second;

      size_t err_count_for_current_index = 0;

      if (!head.hasBox(mapped_box_id)) {
         TBOX_ASSERT(head.hasBox(mapped_box_id));
      }
      const Box& head_mapped_box = *head.getBoxStrict(mapped_box_id);

      for (NeighborSet::OrderedConstIterator na = nabrs.orderedBegin();
           na != nabrs.orderedEnd(); ++na) {

         const Box nabr = *na;

         if (nabr.getOwnerRank() == getBase().getMPI().getRank()) {

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
               for (NeighborSet::OrderedConstIterator nj = nabrs.orderedBegin();
                    nj != nabrs.orderedEnd(); ++nj) {
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
               for (NeighborSet::OrderedConstIterator nj = nabrs.orderedBegin();
                    nj != nabrs.orderedEnd(); ++nj) {
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

            NeighborSet::OrderedConstIterator found_nabr_ =
               nabr_nabrs.find(nabr_nabr);

            if (found_nabr_ == nabr_nabrs.orderedEnd()) {
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
               for (NeighborSet::OrderedConstIterator
                    nj = nabrs.orderedBegin(); nj != nabrs.orderedEnd(); ++nj) {
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
               for (NeighborSet::OrderedConstIterator nj = nabr_nabrs.orderedBegin();
                    nj != nabr_nabrs.orderedEnd(); ++nj) {
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
   left_minus_right.initialize(left.d_base_handle->getBoxLevel(),
      left.d_head_handle->getBoxLevel(),
      left.d_base_width,
      left.getParallelState());
   NeighborhoodSet& drelationships = left_minus_right.d_relationships;

   for (ConstNeighborhoodIterator ai = left.begin(); ai != left.end(); ++ai) {

      const BoxId& mapped_box_id = ai->first;
      const NeighborSet& anabrs = ai->second;

      ConstNeighborhoodIterator bi = right.findLocal(mapped_box_id);
      if (bi != right.end()) {
         const NeighborSet& bnabrs = bi->second;
         // Remove bi from ai.  Put results in a_minus_b.
         NeighborSet& diff = drelationships.getNeighborSet(mapped_box_id,
                                                           left.d_base_width.getDim());

         //equivalent of stl set_difference
         BoxContainer::OrderedConstIterator na = anabrs.orderedBegin();
         BoxContainer::OrderedConstIterator nb = bnabrs.orderedBegin();
         while (na != anabrs.orderedEnd() && nb != bnabrs.orderedEnd()) {
            if ((*na).getId() < (*nb).getId()) {
               diff.insert(diff.orderedEnd(), (*na));
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
         drelationships.clear();
         drelationships.insertNeighborSet(mapped_box_id, anabrs);
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
   const int number_of_inconsistencies =
      checkConsistencyWithHead();
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

   const BoxSet& head_mapped_boxes =
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

         BoxSet::OrderedConstIterator na_in_head =
            head_mapped_boxes.find(unshifted_nabr);

         if (na_in_head == head_mapped_boxes.orderedEnd()) {
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
   BoxList& nbr_boxes) const
{
   const NeighborSet& nbr_mapped_boxes = getNeighborSet(mapped_box_id);
   for (NeighborSet::OrderedConstIterator ni = nbr_mapped_boxes.orderedBegin();
        ni != nbr_mapped_boxes.orderedEnd(); ni++) {
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
   t_initialize_private = tbox::TimerManager::getManager()->
      getTimer("hier::Connector::initializePrivate()");
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
   t_initialize_private.setNull();
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
