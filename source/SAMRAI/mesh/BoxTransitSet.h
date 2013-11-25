/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_BoxTransitSet
#define included_mesh_BoxTransitSet

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/SequentialLocalIdGenerator.h"
#include "SAMRAI/mesh/BalanceBoxBreaker.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <set>

namespace SAMRAI {
namespace mesh {


/*
 * @brief Implementation of TreeLoadBalancer::TransitSet, representing
 * the load with a set of boxes, sorted from highest load to lowest
 * load.
 *
 * As a container, this class is identical to
 * std::set<BoxInTransit,BoxInTransitMoreLoad>.
 */

class BoxTransitSet : public TransitLoad {

public:

   typedef double LoadType;

   BoxTransitSet( const PartitioningParams &pparams );

   //! @brief Initialize implementation for TransitLoad interface.
   void initialize();

   //! @brief Clone object.
   boost::shared_ptr<TransitLoad> clone() const;

   //! @brief Return the total load contained.
   LoadType getSumLoad() const { return d_sumload; }

   //! @brief Insert all boxes from the given BoxContainer.
   void insertAll( const hier::BoxContainer &other );

   //! @brief Insert all boxes from the given BoxTransitSet.
   void insertAll( const TransitLoad &other );

   //! @brief Return number of items in this container.
   size_t getNumberOfItems() const;

   //! @brief Return number of processes contributing to the contents.
   size_t getNumberOfOriginatingProcesses() const;

   //@{
   //! @name Packing/unpacking for communication.

   void putToMessageStream( tbox::MessageStream &msg ) const;

   void getFromMessageStream( tbox::MessageStream &msg );
   //@}


   /*!
    * @brief Adjust the load in this BoxTransitSet by moving work
    * between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoad(
      TransitLoad& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );


   /*!
    * @brief Assign contents to local process and generate
    * balanced<==>unbalanced map.
    *
    * This method uses communication to set up the map.
    */
   void
   assignContentToLocalProcessAndGenerateMap(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector &balanced_to_unbalanced,
      hier::MappingConnector &unbalanced_to_balanced );


   /*
    * @brief Reassign the boxes to the new owner.
    *
    * Any box that isn't already owned by the new owner or doesn't
    * have a valid LocalId, is given one by the
    * SequentialLocalIdGenerator.
    */
   void
   reassignOwnership(
      hier::SequentialLocalIdGenerator &id_gen,
      int new_owner_rank );


   /*!
    * @brief Generate the balanced BoxLevel and most
    * unbalanced<==>balanced edges.  Identify semi-local edges that
    * must be communicated to remote owners.
    */
   void
   generateBalancedBoxLevelAndMostMapEdges(
      hier::BoxLevel &balanced_box_level,
      hier::MappingConnector &unbalanced_to_balanced,
      hier::MappingConnector &balanced_to_unbalanced,
      BoxTransitSet &semi_local ) const;


   //! @brief Whether to allow box breaking.
   void allowBoxBreaking( bool allow_box_breaking ) {
      d_allow_box_breaking = allow_box_breaking;
   }

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named "mesh::BoxTransitSet::*",
    * where the third field is the specific steps performed
    * by the Schedule.  You can override the first two
    * fields with this method.  Conforming to the timer
    * naming convention, timer_prefix should have the form
    * "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   void recursivePrint(
      std::ostream &co=tbox::plog,
      const std::string &border=std::string(),
      int detail_depth=1 ) const;


private:

   static const int BoxTransitSet_EDGETAG0 = 3;
   static const int BoxTransitSet_EDGETAG1 = 4;


   /*!
    * @brief Data to save for each Box that gets passed along the tree
    * edges.
    *
    * The purpose of the BoxInTransit is to associate extra data with
    * a Box as it is broken up and passed from process to process.  A
    * BoxInTransit is a Box going through these changes.  It has a
    * current work load and an orginating Box.
    */
   struct BoxInTransit {

      /*!
       * @brief Constructor
       *
       * @param[in] dim
       */
      BoxInTransit(const tbox::Dimension& dim) :
         d_box(dim),
         d_orig_box(dim) {}


      /*!
       * @brief Construct a new BoxInTransit from an originating box.
       *
       * @param[in] origin
       */
      BoxInTransit( const hier::Box& origin ):
         d_box(origin),
         d_orig_box(origin),
         d_boxload(origin.size()) {}

      /*!
       * @brief Construct new object like an existing object but with a new ID.
       *
       * @param[in] other
       * @param[in] box
       * @param[in] rank
       * @param[in] local_id
       */
      BoxInTransit(
         const BoxInTransit& other,
         const hier::Box& box,
         int rank,
         hier::LocalId local_id) :
         d_box(box, local_id, rank),
         d_orig_box(other.d_orig_box),
         d_boxload(d_box.size()) {}

      /*!
       * @brief Assignment operator
       *
       * @param[in] other
       */
      BoxInTransit& operator = (const BoxInTransit& other) {
         d_box = other.d_box;
         d_orig_box = other.d_orig_box;
         d_boxload = other.d_boxload;
         return *this;
      }

      //! @brief Return the owner rank.
      int getOwnerRank() const {
         return d_box.getOwnerRank();
      }

      //! @brief Return the LocalId.
      hier::LocalId getLocalId() const {
         return d_box.getLocalId();
      }

      //! @brief Return the Box.
      hier::Box& getBox() {
         return d_box;
      }

      //! @brief Return the Box.
      const hier::Box& getBox() const {
         return d_box;
      }

      //! @brief Put self into a MessageStream.
      void putToMessageStream(tbox::MessageStream &mstream) const {
         d_box.putToMessageStream(mstream);
         d_orig_box.putToMessageStream(mstream);
         mstream << d_boxload;
         return;
      }

      //! @brief Set attributes according to data in a MessageStream.
      void getFromMessageStream(tbox::MessageStream &mstream) {
         d_box.getFromMessageStream(mstream);
         d_orig_box.getFromMessageStream(mstream);
         mstream >> d_boxload;
         return;
      }

      hier::Box d_box;

      //! @brief Originating Box (the oldest one leading to this one).
      hier::Box d_orig_box;

      //! @brief Work load in this d_box.
      LoadType d_boxload;
   };


   /*!
    * @brief Comparison functor for sorting BoxInTransit from more to
    * less loads.
    *
    * Ties are broken by BlockId, then lexical comparison of the box's
    * lower corner, then lexical comparison of the upper corner, then
    * orginator BoxId.  The comparison should not use the box's BoxId
    * because some boxes may not have valid ones.
    */
   struct BoxInTransitMoreLoad {
      bool operator () (
         const BoxInTransit& a,
         const BoxInTransit& b) const {
         if (a.getBox().size() != b.getBox().size()) {
            return a.d_boxload > b.d_boxload;
         }
         if ( a.getBox().getBlockId() != b.getBox().getBlockId() ) {
            return a.getBox().getBlockId() < b.getBox().getBlockId();
         }
         if ( a.getBox().lower() != b.getBox().lower() ) {
            return lexicalIndexLessThan(a.getBox().lower(), b.getBox().lower());
         }
         if ( a.getBox().upper() != b.getBox().upper() ) {
            return lexicalIndexLessThan(a.getBox().upper(), b.getBox().upper());
         }
         return a.d_orig_box.getBoxId() < b.d_orig_box.getBoxId();
      }
   private:
      bool lexicalIndexLessThan( const hier::IntVector &a,
                                 const hier::IntVector &b ) const {
         for ( int i=0; i<a.getDim().getValue(); ++i ) {
            if ( a(i) != b(i) ) return a(i) < b(i);
         }
         return false;
      }
   };


   //@{
   //! @name Interfaces like the C++ standard stl::set, to help readability.
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::iterator iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::const_iterator const_iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::reverse_iterator reverse_iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::key_type key_type;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::value_type value_type;
   iterator begin() { return d_set.begin(); }
   iterator end() { return d_set.end(); }
   const_iterator begin() const { return d_set.begin(); }
   const_iterator end() const { return d_set.end(); }
   reverse_iterator rbegin() const { return d_set.rbegin(); }
   reverse_iterator rend() const { return d_set.rend(); }
   size_t size() const { return d_set.size(); }
   std::pair<iterator, bool> insert( const value_type &x ) {
      std::pair<iterator,bool> rval = d_set.insert(x);
      if ( rval.second ) d_sumload += x.d_boxload;
      return rval;
   }
   void erase(iterator pos) { d_sumload -= pos->d_boxload; d_set.erase(pos); }
   size_t erase(const key_type &k) {
      const size_t num_erased = d_set.erase(k);
      if ( num_erased ) d_sumload -= k.d_boxload;
      return num_erased;
   }
   bool empty() const { return d_set.empty(); }
   void clear() { d_sumload = 0; d_set.clear(); }
   void swap( BoxTransitSet &other ) {
      const LoadType tl = d_sumload;
      d_sumload = other.d_sumload;
      other.d_sumload = tl;
      d_set.swap(other.d_set);
   }
   iterator lower_bound( const key_type &k ) const { return d_set.lower_bound(k); }
   iterator upper_bound( const key_type &k ) const { return d_set.upper_bound(k); }
   //@}


   /*!
    * @brief Insert BoxInTransit into an output stream.
    */
   friend std::ostream& operator << ( std::ostream& co,
                                      const BoxInTransit& r) {
      co << r.d_box
         << r.d_box.numberCells() << '|'
         << r.d_box.numberCells().getProduct();
      co << '-' << r.d_orig_box
         << r.d_orig_box.numberCells() << '|'
         << r.d_orig_box.numberCells().getProduct();
      return co;
   }


   //@{ @name Load adjustment methods

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving the
    * biggest between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType adjustLoadByPopping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );

   /*!
    * @brief Adjust the load in this BoxTransitSet by swapping boxes
    * between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType adjustLoadBySwapping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );


   /*!
    * @brief Adjust the load in this BoxTransitSet by moving work
    * between it and another BoxTransitSet.  One box may be broken
    * up to have a part of its load moved.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType adjustLoadByBreaking(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );


   /*!
    * @brief Find a BoxInTransit in each of the source and destination
    * containers that, when swapped, effects a transfer of the given
    * amount of work from the source to the destination.  Swap the boxes.
    *
    * @param [in,out] src
    *
    * @param [in,out] dst
    *
    * @param actual_transfer [out] Amount of work transfered from src to
    * dst.
    *
    * @param ideal_transfer [in] Amount of work to be transfered from
    * src to dst.
    *
    * @param low_transfer
    *
    * @param high_transfer
    */
   bool swapLoadPair(
      BoxTransitSet& src,
      BoxTransitSet& dst,
      LoadType& actual_transfer,
      LoadType ideal_transfer,
      LoadType low_transfer,
      LoadType high_transfer ) const;

   //@}


   /*!
    * @brief Construct semilocal relationships in
    * unbalanced--->balanced Connector.
    *
    * Constructing semilocal unbalanced--->balanced relationships
    * require communication to determine where exported work ended up.
    *
    * @param [out] unbalanced_to_balanced Connector to store
    * relationships in.
    *
    * @param [in] kept_imports Work that was imported and locally kept.
    */
   void constructSemilocalUnbalancedToBalanced(
      hier::MappingConnector &unbalanced_to_balanced,
      const BoxTransitSet &kept_imports ) const;


   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   const BoxTransitSet &recastTransitLoad( const TransitLoad &transit_load ) {
      const BoxTransitSet *ptr = static_cast<const BoxTransitSet*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   BoxTransitSet &recastTransitLoad( TransitLoad &transit_load ) {
      BoxTransitSet *ptr = static_cast<BoxTransitSet*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void initializeCallback() {
      TimerStruct& timers(s_static_timers[s_default_timer_prefix]);
      getAllTimers(s_default_timer_prefix, timers);
   }


   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void finalizeCallback() {
      s_static_timers.clear();
   }


   //! @brief Compute the load for a Box.
   double computeLoad( const hier::Box& box) const {
      return double(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      return double((box * restriction).size());
   }

   /*!
    * @brief Look for an input database called "BoxTransitSet" and
    * read parameters if it exists.
    */
   void
   getFromInput();


   //! @brief Balance penalty is proportional to imbalance.
   double computeBalancePenalty( double imbalance) const {
      return tbox::MathUtilities<double>::Abs(imbalance);
   }


   std::set<BoxInTransit, BoxInTransitMoreLoad> d_set;
   LoadType d_sumload;

   const PartitioningParams *d_pparams;

   BalanceBoxBreaker d_box_breaker;

   bool d_allow_box_breaking;


   //@{
   //! @name Debugging stuff.
   bool d_print_steps;
   bool d_print_pop_steps;
   bool d_print_swap_steps;
   bool d_print_break_steps;
   bool d_print_edge_steps;
   //@}

   //@{
   //! @name Timer data for Schedule class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each Schedule object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      boost::shared_ptr<tbox::Timer> t_adjust_load;
      boost::shared_ptr<tbox::Timer> t_adjust_load_by_popping;
      boost::shared_ptr<tbox::Timer> t_adjust_load_by_swapping;
      boost::shared_ptr<tbox::Timer> t_shift_loads_by_breaking;
      boost::shared_ptr<tbox::Timer> t_find_swap_pair;
      boost::shared_ptr<tbox::Timer> t_assign_content_to_local_process_and_generate_map;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal_comm_wait;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal_send_edges;
      boost::shared_ptr<tbox::Timer> t_pack_edge;
      boost::shared_ptr<tbox::Timer> t_unpack_edge;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   /*!
    * @brief Structure of timers in s_static_timers, matching this
    * object's timer prefix.
    */
   TimerStruct* d_object_timers;

   /*!
    * @brief Get all the timers defined in TimerStruct.  The timers
    * are named with the given prefix.
    */
   static void
   getAllTimers(
      const std::string& timer_prefix,
      TimerStruct& timers);

   //@}

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};



}
}

#endif
