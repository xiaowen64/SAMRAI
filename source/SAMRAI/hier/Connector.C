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
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealMappedBoxConstIterator.h"
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

tbox::Pointer<tbox::Timer> Connector::t_initialize;

tbox::Pointer<tbox::Timer> Connector::t_acquire_remote_relationships;

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
   d_relationships(),
   d_global_relationships(),
   d_parallel_state(MappedBoxLevel::DISTRIBUTED),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_rank(BAD_INT),
   d_nproc(BAD_INT),
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
   d_parallel_state(other.d_parallel_state),
   d_global_number_of_neighbor_sets(other.d_global_number_of_neighbor_sets),
   d_global_number_of_relationships(other.d_global_number_of_relationships),
   d_global_data_up_to_date(other.d_global_data_up_to_date),
   d_rank(other.d_rank),
   d_nproc(other.d_nproc),
   d_connector_type(other.d_connector_type)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const MappedBoxLevel& base_mapped_box_level,
   const MappedBoxLevel& head_mapped_box_level,
   const IntVector& base_width,
   const NeighborhoodSet& relationships,
   const MappedBoxLevel::ParallelState parallel_state)
   :d_base_width(base_width.getDim(), 0),
   d_ratio(base_width.getDim(), 0),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_parallel_state(MappedBoxLevel::DISTRIBUTED),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_rank(BAD_INT),
   d_nproc(BAD_INT),
   d_connector_type(UNKNOWN)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(base_mapped_box_level,
      head_mapped_box_level,
      base_width);

   initialize(base_mapped_box_level,
      head_mapped_box_level,
      base_width,
      relationships,
      parallel_state);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const MappedBoxLevel& base_mapped_box_level,
   const MappedBoxLevel& head_mapped_box_level,
   const IntVector& base_width,
   const MappedBoxLevel::ParallelState parallel_state)
   :d_base_width(base_width.getDim(), 0),
   d_ratio(base_width.getDim(), 0),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_parallel_state(MappedBoxLevel::DISTRIBUTED),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(true),
   d_rank(BAD_INT),
   d_nproc(BAD_INT),
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
& Connector::getNeighborhoodSets() const
{
   TBOX_ASSERT(isInitialized());
   return d_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

const NeighborhoodSet
& Connector::getGlobalNeighborhoodSets() const
{
   TBOX_ASSERT(isInitialized());
   if (d_parallel_state == MappedBoxLevel::DISTRIBUTED) {
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
   const MappedBoxId& mapped_box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == MappedBoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != d_rank) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasMappedBox(mapped_box_id)) {
      TBOX_ERROR(
         "Exiting due to above reported error."
         << "Connector::insertNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == MappedBoxLevel::GLOBALIZED) {
      d_global_relationships[mapped_box_id].insert(neighbors.begin(),
         neighbors.end());
   }
   if (mapped_box_id.getOwnerRank() == d_rank) {
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
   const MappedBoxId& mapped_box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == MappedBoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != d_rank) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasMappedBox(mapped_box_id)) {
      TBOX_ERROR(
         "Connector::eraseNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == MappedBoxLevel::GLOBALIZED) {
      NeighborhoodSet::iterator mi(d_global_relationships.find(mapped_box_id));
      if (mi != d_global_relationships.end()) {
         mi->second.erase(neighbor);
      }
   }
   if (mapped_box_id.getOwnerRank() == d_rank) {
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
void Connector::swapNeighbors(
   NeighborSet& neighbors,
   const MappedBoxId& mapped_box_id)
{
   TBOX_ASSERT(mapped_box_id.getPeriodicId() == PeriodicId::zero());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == MappedBoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != d_rank) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasMappedBox(mapped_box_id)) {
      TBOX_ERROR(
         "Exiting due to above reported error."
         << "Connector::swapNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (mapped_box_id.getOwnerRank() == d_rank) {
      if (d_parallel_state == MappedBoxLevel::GLOBALIZED) {
         d_global_relationships[mapped_box_id] = neighbors;
      }
      d_relationships[mapped_box_id].swap(neighbors);
   } else if (d_parallel_state == MappedBoxLevel::GLOBALIZED) {
      d_global_relationships[mapped_box_id].swap(neighbors);
   }
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

   std::vector<int> recv_mesg_size(d_nproc);
   mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(d_nproc);
   int totl_size = 0;
   for (int n = 0; n < d_nproc; ++n) {
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
   const int num_mapped_boxes = static_cast<int>(d_relationships.size());
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

      const MappedBoxId& mapped_box_id = (*ci).first;
      const NeighborSet& nabrs = (*ci).second;

      send_mesg[imesg++] = mapped_box_id.getLocalId().getValue();
      send_mesg[imesg++] = mapped_box_id.getBlockId().getBlockValue();
      send_mesg[imesg++] = static_cast<int>(nabrs.size());

      for (NeighborSet::const_iterator ni = nabrs.begin();
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

   for (int n = 0; n < d_nproc; ++n) {

      if (n != d_rank) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int num_nabr_lists = *(ptr++);

         for (int i = 0; i < num_nabr_lists; ++i) {

            const LocalId local_id(*ptr++);
            const BlockId block_id(*ptr++);
            const int num_nabrs = (*ptr++);
            const MappedBoxId mapped_box_id(local_id, n, block_id);

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
   const MappedBoxLevel::ParallelState parallel_state)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isInitialized()) {
      TBOX_ERROR(
         "Connector::setParallelState: Cannot change the parallel state of\n"
         << "an uninitialized Connector.  See Connector::initialize()");
   }
#endif
   if (parallel_state != MappedBoxLevel::DISTRIBUTED && parallel_state != MappedBoxLevel::GLOBALIZED) {
      TBOX_ERROR("Connector::setParallelState: Invalid distribution state: "
         << parallel_state << "\n");
   }

   if (d_parallel_state == MappedBoxLevel::DISTRIBUTED && parallel_state == MappedBoxLevel::GLOBALIZED) {
      acquireRemoteNeighborhoods();
   } else if (d_parallel_state == MappedBoxLevel::GLOBALIZED && parallel_state ==
              MappedBoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initialize(
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_width,
   const NeighborhoodSet& relationships,
   const MappedBoxLevel::ParallelState parallel_state)
{
   if (&relationships != &d_relationships) {
      d_relationships = relationships;
   }
   initializePrivate(base,
      head,
      base_width,
      base.getMappedBoxLevelHandle()->getMappedBoxLevel().getRefinementRatio(),
      head.getMappedBoxLevelHandle()->getMappedBoxLevel().getRefinementRatio(),
      parallel_state);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initialize(
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_width,
   const MappedBoxLevel::ParallelState parallel_state)
{
   NeighborhoodSet dummy_relationships;
   swapInitialize(base, head, base_width, dummy_relationships, parallel_state);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::swapInitialize(
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_width,
   NeighborhoodSet& relationships,
   const MappedBoxLevel::ParallelState parallel_state)
{
   TBOX_ASSERT(&relationships != &d_relationships);   // Library error if this fails.
   d_relationships.swap(relationships);
   initializePrivate(base,
      head,
      base_width,
      base.getRefinementRatio(),
      head.getRefinementRatio(),
      parallel_state);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::initializePrivate(
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_width,
   const IntVector& baseRefinementRatio,
   const IntVector& headRefinementRatio,
   const MappedBoxLevel::ParallelState parallel_state)
{
   t_initialize->start();
   /*
    * Check inputs.
    */
   if (!base.isInitialized() ||
       !head.isInitialized()) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Connector may not be initialized with\n"
         << "an uninitialized MappedBoxLevel.");
   }
   if (base.getGridGeometry() != head.getGridGeometry()) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Connector must be initialized with\n"
         << "MappedBoxLevels using the same GridGeometry.");
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
   if (!(base_width >= hier::IntVector::getZero(base_width.getDim()))) {
      TBOX_ERROR("Connector::initializePrivate():\n"
         << "Invalid ghost cell width: "
         << base_width << "\n");
   }
   if (parallel_state == MappedBoxLevel::GLOBALIZED &&
       base.getParallelState() != MappedBoxLevel::GLOBALIZED) {
      TBOX_ERROR(
         "Connector::initializePrivate: base MappedBoxLevel must be in\n"
         << "GLOBALIZED state before initializing the Connector to\n"
         << "GLOBALIZED state.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   bool errf = false;
   for (NeighborhoodSet::const_iterator ci = d_relationships.begin(); ci != d_relationships.end(); ++ci) {
      if (!base.hasMappedBox((*ci).first)) {
         const NeighborSet& nabrs = (*ci).second;
         tbox::perr << "\nConnector::initializePrivate: NeighborhoodSet "
                    << "provided for non-existent mapped_box " << ci->first
                    << "\n" << "Neighbors (" << nabrs.size() << "):\n";
         for (NeighborSet::const_iterator na = nabrs.begin();
              na != nabrs.end(); ++na) {
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

   d_nproc = base.getNproc();
   d_rank = base.getRank();

   d_base_handle = base.getMappedBoxLevelHandle();
   d_head_handle = head.getMappedBoxLevelHandle();

   d_base_width = base_width;

   if (baseRefinementRatio <= headRefinementRatio) {
      d_ratio = headRefinementRatio / baseRefinementRatio;
      d_head_coarser = false;
      d_ratio_is_exact = ( d_ratio * baseRefinementRatio ) == headRefinementRatio;
   } else {
      d_ratio = baseRefinementRatio / headRefinementRatio;
      d_head_coarser = true;
      d_ratio_is_exact = ( d_ratio * headRefinementRatio ) == baseRefinementRatio;
   }
   if (baseRefinementRatio * headRefinementRatio <
       IntVector::getZero(base_width.getDim())) {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      d_ratio = -headRefinementRatio * baseRefinementRatio;
      d_ratio_is_exact = true;
   }

   if (parallel_state == MappedBoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   } else {
      if (&d_relationships != &d_global_relationships) {
         d_global_relationships = d_relationships;
      }
   }

   // Erase remote relationships, if any, from d_relationships.
   if (!d_relationships.empty()) {
      if (d_relationships.begin()->first.getOwnerRank() != d_rank ||
          d_relationships.rbegin()->first.getOwnerRank() != d_rank) {
         NeighborhoodSet::Range range = d_relationships.findRanksRange(d_rank);
         if (range.first == range.second) {
            // No relationship belongs to local process.
            d_relationships.clear();
         } else {
            // Erase relationships belonging to remote process < d_rank.
            d_relationships.erase(d_relationships.begin(), range.first);
            // Erase relationships belonging to remote process > d_rank.
            d_relationships.erase(range.second, d_relationships.end());
         }
      }
   }

   d_parallel_state = parallel_state;

   t_initialize->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::clear()
{
   if (d_nproc != BAD_INT) {
      d_relationships.clear();
      d_global_relationships.clear();
      d_base_handle.setNull();
      d_head_handle.setNull();
      d_base_width(0) = d_ratio(0) = 0;
      d_parallel_state = MappedBoxLevel::DISTRIBUTED;
      d_rank = d_nproc = BAD_INT;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void Connector::swap(
   Connector& a,
   Connector& b)
{

   if (&a != &b) {
      tbox::Pointer<MappedBoxLevelHandle> tmplayer;
      // tbox::Pointer<PersistentOverlapConnectors> tmphub;
      int tmpint;
      bool tmpbool;
      IntVector tmpvec(a.getBase().getDim());
      MappedBoxLevel::ParallelState tmpstate;

      tmplayer = a.d_base_handle;
      a.d_base_handle = b.d_base_handle;
      b.d_base_handle = tmplayer;

      tmplayer = a.d_head_handle;
      a.d_head_handle = b.d_head_handle;
      b.d_head_handle = tmplayer;

      a.d_relationships.swap(b.d_relationships);
      a.d_global_relationships.swap(b.d_global_relationships);

      tmpvec = a.d_ratio;
      a.d_ratio = b.d_ratio;
      b.d_ratio = tmpvec;

      tmpbool = a.d_head_coarser;
      a.d_head_coarser = b.d_head_coarser;
      b.d_head_coarser = tmpbool;

      tmpstate = a.d_parallel_state;
      a.d_parallel_state = b.d_parallel_state;
      b.d_parallel_state = tmpstate;

      tmpint = a.d_nproc;
      a.d_nproc = b.d_nproc;
      b.d_nproc = tmpint;

      tmpint = a.d_rank;
      a.d_rank = b.d_rank;
      b.d_rank = tmpint;

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
      connector.d_head_handle->getMappedBoxLevel(),
      connector.d_base_handle->getMappedBoxLevel(),
      my_gcw,
      MappedBoxLevel::DISTRIBUTED);
   TBOX_ASSERT(isTransposeOf(connector));

   const tbox::Dimension dim(my_gcw.getDim());
   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   const NeighborhoodSet& r_relationships = connector.getNeighborhoodSets();

   for (NeighborhoodSet::const_iterator ci = r_relationships.begin();
        ci != r_relationships.end(); ++ci) {

      const MappedBoxId& mapped_box_id = ci->first;
      const MappedBoxSet::const_iterator ni = getHead().getMappedBox(mapped_box_id);
      if (ni == getHead().getMappedBoxes().end()) {
         TBOX_ERROR(
            "Connector::initializeToLocalTranspose: mapped_box index\n"
            << mapped_box_id
            << " not found in local part of head mapped_box_level.\n"
            << "This means that the incoming Connector data was not a\n"
            << "self-consistent local mapping.\n");
      }
      const Box& my_head_mapped_box = *ni;

      const NeighborSet& my_base_subset = ci->second;
      for (NeighborSet::const_iterator na = my_base_subset.begin();
           na != my_base_subset.end(); ++na) {
         const Box& my_base_mapped_box = *na;
         if (my_base_mapped_box.getOwnerRank() != d_rank) {
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
            if (getHead().hasMappedBox(my_shifted_head_mapped_box)) {
               hier::MappedBoxId base_non_per_id(
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
      if ((*ei).second.empty()) {
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

bool
Connector::isTransposeOf(
   const Connector& other) const
{
   bool rval = false;
   if (d_base_handle == other.d_head_handle && d_head_handle == other.d_base_handle) {
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
      const MappedBoxSet& nabrs = ei->second;
      for (MappedBoxSet::const_iterator na = nabrs.begin();
           na != nabrs.end();
           ++na) {
         if (na->getOwnerRank() != d_rank) {
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
   // FIXME:  BUG
   TBOX_ASSERT(isInitialized());
   size_t local_number_of_relationships = 0;
   for (NeighborhoodSet::const_iterator ei(d_relationships.begin()); ei != d_relationships.end(); ++ei) {
      d_global_number_of_relationships += static_cast<int>(ei->second.size());
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
   tbox::SAMRAI_MPI mpi(getMPI());
   TBOX_ASSERT(isInitialized());

   if (d_global_data_up_to_date) {
      return;
   }
   if (d_parallel_state == MappedBoxLevel::GLOBALIZED) {
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
      (head_gcw * ratio) : IntVector::ceiling(head_gcw, ratio);

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
   const NeighborhoodSet& relationships = getNeighborhoodSets();
   const hier::IntVector head_gcw =
      convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         d_base_width);
   os << border << "Parallel state     : "
      << (getParallelState() == MappedBoxLevel::DISTRIBUTED ? "DIST" : "GLOB")
      << '\n'
      << border << "Rank,nproc         : " << d_rank << ", " << d_nproc << '\n'
      << border << "Base,head objects  :"
      << " ("
      << (d_base_handle == d_head_handle ? "same" : "different") << ") "
      << (void *)&d_base_handle->getMappedBoxLevel() << ", "
      << (void *)&d_head_handle->getMappedBoxLevel() << "\n"
      << border << "Base,head,/ ratios : "
      << getBase().getRefinementRatio() << ", "
      << getHead().getRefinementRatio() << ", "
      << d_ratio << (d_head_coarser ? " (head coarser)" : "") << '\n'
      << border << "Base,head widths   : " << d_base_width << ", "
      << head_gcw << '\n'
      << border << "Box count    : " << getBase().getLocalNumberOfBoxes()
      << " (" << relationships.size() << " with neighbor lists)\n"
   ;
   if (detail_depth > 0) {
      os << border << "Mapped_boxes with neighbors:\n";
      for (NeighborhoodSet::const_iterator ei = relationships.begin(); ei != relationships.end();
           ++ei) {
         MappedBoxSet::const_iterator ni = getBase().getMappedBox(ei->first);
         if (ni != getBase().getMappedBoxes().end()) {
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
         hier::Box ghost_box = (*ni);
         ghost_box.grow(d_base_width);
         const NeighborSet& nabrs = (*ei).second;
         os << border << "    Neighbors (" << nabrs.size() << "):"
            << ((detail_depth > 1) ? "\n" : " ...\n");
         if (detail_depth > 1) {
            NeighborSet::const_iterator i_nabr;
            for (i_nabr = nabrs.begin(); i_nabr != nabrs.end(); ++i_nabr) {
               hier::Box ovlap = *i_nabr;
               if ( ni->getBlockId() != i_nabr->getBlockId() ) {
                  d_base_handle->getMappedBoxLevel().getGridGeometry()->
                     translateBox(
                        ovlap,
                        d_head_handle->getMappedBoxLevel().getRefinementRatio(),
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

   const NeighborhoodSet& relationships = getNeighborhoodSets();

   for (NeighborhoodSet::const_iterator ni = relationships.begin();
        ni != relationships.end(); ++ni) {

      ++sum_nabr_sets;
      const Box& mapped_box = *getBase().getMappedBoxStrict(ni->first);
      sum_cells += mapped_box.size();
      Box base_ghost_box = mapped_box;
      base_ghost_box.grow(d_base_width);

      const NeighborSet& nabrs = ni->second;

      max_nabrs = tbox::MathUtilities<int>::Max(max_nabrs, (int)nabrs.size());
      min_nabrs = tbox::MathUtilities<int>::Min(min_nabrs, (int)nabrs.size());

      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {

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

         if (nabr.getOwnerRank() == d_rank) {
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
   if (d_nproc > 1) {
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
      << ' ' << std::setw(8) << std::right << static_cast<float>(sum_num[i] / d_nproc)
      << '\n';
   }

   // Write statistics that cannot be expressed as average over P processors or N MappedBoxes.
   os << border
      <<
   "                               local        min               max             avg\n";
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
   const Connector &connector,
   const std::string& border,
   int detail_depth )
   : d_conn(connector),
     d_border(border),
     d_detail_depth(detail_depth)
{
   return;
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
   format.d_conn.recursivePrint( os, format.d_border, format.d_detail_depth );
   return os;
}


/*
 ***********************************************************************
 * Return a Outputter that can dump the Connector to a stream.
 ***********************************************************************
 */

Connector::Outputter Connector::format(
   const std::string& border,
   int detail_depth ) const
{
   return Outputter( *this, border, detail_depth);
}


/*
 ***********************************************************************
 ***********************************************************************
 */

Connector *Connector::makeGlobalizedCopy(
   const Connector& other) const
{
   // Prevent wasteful accidental use when this method is not needed.
   TBOX_ASSERT(other.getParallelState() != MappedBoxLevel::GLOBALIZED);

   Connector* copy = new Connector(other);
   copy->setParallelState(MappedBoxLevel::GLOBALIZED);
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
      (input_transpose.d_parallel_state == MappedBoxLevel::GLOBALIZED) ?
      &input_transpose : makeGlobalizedCopy(input_transpose);

   const MappedBoxLevel &head = getHead().getGlobalizedVersion();

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   /*
    * Check for extraneous relationships.
    * For every relationship in this, there should be reverse relationship in transpose.
    */
   const NeighborhoodSet& this_relationships = getNeighborhoodSets();

   Box shifted_mapped_box(dim);   // Shifted version of an unshifted Box.
   Box unshifted_mapped_box(dim); // Unhifted version of a shifted Box.

   size_t err_count = 0;

   for (NeighborhoodSet::const_iterator ci = this_relationships.begin();
        ci != this_relationships.end(); ++ci) {

      const MappedBoxId& mapped_box_id = ci->first;
      const Box& mapped_box = *getBase().getMappedBox(mapped_box_id);

      size_t err_count_for_current_index = 0;

      const NeighborSet& nabrs = ci->second;

      for (NeighborSet::const_iterator ni = nabrs.begin();
           ni != nabrs.end(); ++ni) {

         if (ignore_periodic_relationships && ni->isPeriodicImage()) {
            continue;
         }

         const Box& nabr = *ni;

         const NeighborhoodSet& tran_relationships = transpose->getGlobalNeighborhoodSets();

         /*
          * Key for find in NeighborhoodSet must be non-periodic.
          */
         MappedBoxId non_per_nabr_id(nabr.getGlobalId(),
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

         NeighborSet::const_iterator nabr_ni;

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
            for (NeighborSet::const_iterator nj = nabr_nabrs.begin();
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
         for (NeighborSet::const_iterator nj = nabrs.begin();
              nj != nabrs.end(); ++nj) {
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

      const MappedBoxId& mapped_box_id = ci->first;
      const NeighborSet& nabrs = ci->second;

      size_t err_count_for_current_index = 0;

      if (!head.hasMappedBox(mapped_box_id)) {
         TBOX_ASSERT(head.hasMappedBox(mapped_box_id));
      }
      const Box& head_mapped_box = *head.getMappedBoxStrict(mapped_box_id);

      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {

         const Box nabr = *na;

         if (nabr.getOwnerRank() == d_rank) {

            if (ignore_periodic_relationships && nabr.isPeriodicImage()) {
               continue;
            }

            if (!getBase().hasMappedBox(nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
                          << "Head mapped_box " << head_mapped_box
                          << " has neighbor " << nabr << "\n"
                          << " but the neighbor does not exist "
                          << "in the base mapped_box_level.\n";
               tbox::perr << "Neighbors of head mapped_box "
                          << mapped_box_id << " are:\n";
               for (NeighborSet::const_iterator nj = nabrs.begin();
                    nj != nabrs.end(); ++nj) {
                  tbox::perr << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box& base_mapped_box = *getBase().getMappedBoxStrict(
                  nabr);

            /*
             * Non-periodic MappedBoxId needed for NeighborhoodSet::find()
             */
            MappedBoxId base_non_per_id(base_mapped_box.getGlobalId(),
                                        base_mapped_box.getBlockId(),
                                        PeriodicId::zero());
            NeighborhoodSet::const_iterator nabr_nabrs_ =
               this_relationships.find(base_non_per_id);

            if (nabr_nabrs_ == this_relationships.end()) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
                          << "Head mapped_box " << head_mapped_box << "\n"
                          << " has base mapped_box "
                          << base_mapped_box << " as a neighbor.\n"
                          << "But " << base_mapped_box
                          << " has no neighbor container.\n";
               tbox::perr << "Neighbors of head mapped_box " << MappedBoxId(
                  mapped_box_id)
                          << ":" << std::endl;
               for (NeighborSet::const_iterator nj = nabrs.begin();
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

            NeighborSet::const_iterator found_nabr_ =
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
               for (NeighborSet::const_iterator
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
               for (NeighborSet::const_iterator nj = nabr_nabrs.begin();
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
   const NeighborhoodSet& relationships = getNeighborhoodSets();
   NeighborhoodSet::const_iterator i_relationships;
   for (i_relationships = relationships.begin(); i_relationships != relationships.end(); ++i_relationships) {
      const MappedBoxId& mapped_box_id = (*i_relationships).first;
      if (!getBase().hasMappedBox(mapped_box_id)) {
         ++num_errors;
         tbox::plog << "ERROR->"
                    <<
         "Connector::assertConsistencyWithBase: Neighbor data given "
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
   if (0) {
      tbox::plog << "Computing relationship differences, a:\n" << left.format(dbgbord, 3)
                 << "Computing relationship differences, b:\n" << right.format(dbgbord, 3);
   }
   left_minus_right.initialize(left.d_base_handle->getMappedBoxLevel(),
      left.d_head_handle->getMappedBoxLevel(),
      left.d_base_width,
      left.getParallelState());
   const NeighborhoodSet& arelationships = left.getNeighborhoodSets();
   const NeighborhoodSet& brelationships = right.getNeighborhoodSets();
   NeighborhoodSet& drelationships = left_minus_right.d_relationships;

   for (NeighborhoodSet::const_iterator ai = arelationships.begin(); ai != arelationships.end();
        ++ai) {

      const MappedBoxId& mapped_box_id = ai->first;
      const NeighborSet& anabrs = ai->second;

      NeighborhoodSet::const_iterator bi = brelationships.find(mapped_box_id);
      if (bi != brelationships.end()) {
         const NeighborSet& bnabrs = (*bi).second;
         // Remove bi from ai.  Put results in a_minus_b.
         NeighborSet& diff = drelationships[mapped_box_id];
         std::insert_iterator<NeighborSet> ii(diff, diff.begin());
         set_difference(anabrs.begin(),
            anabrs.end(),
            bnabrs.begin(),
            bnabrs.end(),
            ii, Box::id_less());
         if (diff.empty()) {
            drelationships.erase(mapped_box_id);
         }
      } else if (!anabrs.empty()) {
         drelationships[mapped_box_id] = anabrs;
      }

   }
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
   const MappedBoxLevel &globalized_head = getHead().getGlobalizedVersion();

   const int number_of_inconsistencies =
      static_cast<int>(
         checkConsistencyWithHead(
            getNeighborhoodSets(),
            globalized_head));

   return number_of_inconsistencies;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t Connector::checkConsistencyWithHead(
   const NeighborhoodSet& relationships,
   const MappedBoxLevel& head_mapped_box_level)
{
   TBOX_ASSERT(head_mapped_box_level.getParallelState() ==
      MappedBoxLevel::GLOBALIZED);

   const MappedBoxSet& head_mapped_boxes =
      head_mapped_box_level.getGlobalMappedBoxes();

   size_t number_of_inconsistencies = 0;

   /*
    * For each neighbor in each neighbor list,
    * check that the neighbor is in the head_mapped_box_level.
    */

   for (NeighborhoodSet::const_iterator ei = relationships.begin();
        ei != relationships.end(); ++ei) {

      const MappedBoxId& mapped_box_id = ei->first;
      const NeighborSet& nabrs = ei->second;

      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {

         const Box& nabr = *na;
         const Box unshifted_nabr(
            nabr, PeriodicId::zero(), head_mapped_box_level.getRefinementRatio());

         MappedBoxSet::const_iterator na_in_head =
            head_mapped_boxes.find(unshifted_nabr);

         if (na_in_head == head_mapped_boxes.end()) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
                       << "Neighbor list for mapped_box " << mapped_box_id << "\n"
                       << "referenced nonexistent neighbor "
                       << nabr << "\n";
            tbox::perr << "Neighbors of mapped_box " << mapped_box_id << ":\n";
            for (MappedBoxSet::const_iterator nb = nabrs.begin();
                 nb != nabrs.end(); ++nb) {
               tbox::perr << "    " << *nb << '\n';
            }
            ++number_of_inconsistencies;
            continue;
         }

         const Box& nabr_in_head = *na_in_head;
         if (! unshifted_nabr.isIdEqual(nabr_in_head) ||
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
   const MappedBoxId& mapped_box_id,
   BoxList& nbr_boxes) const
{
   const NeighborSet& nbr_mapped_boxes = getNeighborSet(mapped_box_id);
   for (NeighborSet::const_iterator ni = nbr_mapped_boxes.begin();
        ni != nbr_mapped_boxes.end(); ni++) {
      nbr_boxes.appendItem(*ni);
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
   t_initialize = tbox::TimerManager::getManager()->
      getTimer("hier::Connector::initialize()");
   t_acquire_remote_relationships = tbox::TimerManager::getManager()->
      getTimer("hier::Connector::acquire_remote_relationships");
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
   t_initialize.setNull();
   t_acquire_remote_relationships.setNull();
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
