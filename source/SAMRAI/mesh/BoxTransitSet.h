/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer::TransitSet.
 *
 ************************************************************************/

#ifndef included_mesh_BoxTransitSet
#define included_mesh_BoxTransitSet

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/SequentialLocalIdGenerator.h"
#include "SAMRAI/mesh/BalanceBoxBreaker.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MessageStream.h"

#include <set>

namespace SAMRAI {
namespace mesh {


/*
 * @brief Implementation of TreeLoadBalancer::TransitSet, representing
 * the load with a set of boxes.
 */

class BoxTransitSet {

public:

   typedef double LoadType;

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
         d_orig_box(dim)
         { }


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
       * @brief Construct new object having the history an existing
       * object but is otherwise different.
       *
       * @param[in] other
       *
       * @param[in] box
       *
       * @param[in] rank
       *
       * @param[in] local_id
       */
      BoxInTransit(
         const BoxInTransit& other,
         const hier::Box& box,
         int rank,
         hier::LocalId local_id) :
         d_box(box, local_id, rank),
         d_orig_box(other.d_orig_box),
         d_boxload(d_box.size())
         { }

      /*!
       * @brief Assignment operator
       *
       * @param[in] other
       */
      const BoxInTransit&
      operator = (const BoxInTransit& other)
      {
         d_box = other.d_box;
         d_orig_box = other.d_orig_box;
         d_boxload = other.d_boxload;
         return *this;
      }

      //! @brief Return the owner rank.
      int
      getOwnerRank() const
      {
         return d_box.getOwnerRank();
      }

      //! @brief Return the LocalId.
      hier::LocalId
      getLocalId() const
      {
         return d_box.getLocalId();
      }

      //! @brief Return the Box.
      hier::Box&
      getBox()
      {
         return d_box;
      }

      //! @brief Return the Box.
      const hier::Box&
      getBox() const
      {
         return d_box;
      }

      /*!
       * @brief Put self into a MessageStream.
       *
       * This is the opposite of getFromMessageStream().
       */
      void
      putToMessageStream(tbox::MessageStream &mstream) const
         {
            d_box.putToMessageStream(mstream);
            d_orig_box.putToMessageStream(mstream);
            mstream << d_boxload;
            return;
         }

      /*!
       * @brief Set attributes according to data in a MessageStream.
       *
       * This is the opposite of putToMessageStream().
       */
      void
      getFromMessageStream(tbox::MessageStream &mstream)
         {
            d_box.getFromMessageStream(mstream);
            d_orig_box.getFromMessageStream(mstream);
            mstream >> d_boxload;
            return;
         }

      //! @brief The Box.
      hier::Box d_box;

      //! @brief Originating Box.
      hier::Box d_orig_box;

      //! @brief Work load in this box.
      LoadType d_boxload;
   };


   /*!
    * @brief Insert BoxInTransit into an output stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BoxInTransit& r)
      {
         co << r.d_box
            << r.d_box.numberCells() << '|'
            << r.d_box.numberCells().getProduct();
         co << '-' << r.d_orig_box
            << r.d_box.numberCells() << '|'
            << r.d_box.numberCells().getProduct();
         return co;
      }


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
      bool
      operator () (
         const BoxInTransit& a,
         const BoxInTransit& b) const
      {
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
      bool lexicalIndexLessThan( const hier::IntVector &a,
                                 const hier::IntVector &b ) const {
         for ( int i=0; i<a.getDim().getValue(); ++i ) {
            if ( a(i) != b(i) ) return a(i) < b(i);
         }
         return false;
      }
   };



public:

   //@{ @name Load adjustment

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving work
    * between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin
    *
    * @param[in] ideal_load The load that main_bin should have.
    *
    * @param[in] low_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoad(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving the
    * biggest between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin
    *
    * @param[in] ideal_load The load that main_bin should have.
    *
    * @param[in] low_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadByPopping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );

   /*!
    * @brief Adjust the load in this BoxTransitSet by swapping boxes
    * between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin
    *
    * @param[in] ideal_load The load that main_bin should have.
    *
    * @param[in] low_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadBySwapping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving work
    * between it and another BoxTransitSet.  One box may be broken
    * up to have a part of its load moved.
    *
    * @param[in,out] hold_bin
    *
    * @param[in] ideal_load The load that main_bin should have.
    *
    * @param[in] low_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when main_bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadByBreaking(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );

   /*!
    * @brief Allow box breaking when adjusting load.
    */
   void allowBoxBreaking() {
      d_allow_box_breaking = true;
   }

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
   bool
   swapLoadPair(
      BoxTransitSet& src,
      BoxTransitSet& dst,
      LoadType& actual_transfer,
      LoadType ideal_transfer,
      LoadType low_transfer,
      LoadType high_transfer ) const;

   //@}

   /*!
    * @brief A set of BoxInTransit, sorted from highest load to lowest load.
    *
    * This class is identical to std::set<BoxInTransit,BoxInTransitMoreLoad>
    * and adds tracking of the sum of loads in the set.
    */
   public:
      //@{
      //! @name Duplicated set interfaces.
      typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::iterator iterator;
      typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::const_iterator const_iterator;
      typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::reverse_iterator reverse_iterator;
      typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::key_type key_type;
      typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::value_type value_type;
   BoxTransitSet();

   void setPartitioningParams( const PartitioningParams &pparams )
      {
         d_pparams = &pparams;
         d_bbb.setPartitioningParams(pparams);
      }
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
      LoadType getSumLoad() const { return d_sumload; }
      void insertAll( const hier::BoxContainer &other ) {
         size_t old_size = d_set.size();
         for ( hier::BoxContainer::const_iterator bi=other.begin(); bi!=other.end(); ++bi ) {
            BoxInTransit new_box(*bi);
            d_set.insert( new_box );
            d_sumload += new_box.d_boxload;
         };
         if ( d_set.size() != old_size + other.size() ) {
            TBOX_ERROR("BoxTransitSet's insertAll currently can't weed out duplicates.");
         }
      }
      void insertAll( const BoxTransitSet &other ) {
         size_t old_size = d_set.size();
         d_set.insert( other.d_set.begin(), other.d_set.end() );
         d_sumload += other.d_sumload;
         if ( d_set.size() != old_size + other.size() ) {
            TBOX_ERROR("BoxTransitSet's insertAll currently can't weed out duplicates.");
         }
      }
      size_t getNumberOfOriginatingProcesses() const {
         std::set<int> originating_procs;
         for ( const_iterator si=begin(); si!=end(); ++si ) {
            originating_procs.insert( si->d_orig_box.getOwnerRank() );
         }
         return originating_procs.size();
      }

   //@{
   //! @name Packing/unpacking for communication.
   void putToMessageStream( tbox::MessageStream &msg ) const;
   void getFromMessageStream( tbox::MessageStream &msg );
   //@}

   private:

   void setTimers();

   /*!
    * @brief Compute the load for a Box.
    */
   double
   computeLoad(
      const hier::Box& box) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the cell.
       */
      return double(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double
   computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the overlap region.
       */
      return double((box * restriction).size());
   }

   double
   computeBalancePenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b,
      double imbalance) const
   {
      NULL_USE(a);
      NULL_USE(b);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   double
   computeBalancePenalty(
      const BoxTransitSet& a,
      const BoxTransitSet& b,
      double imbalance) const
   {
      NULL_USE(a);
      NULL_USE(b);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   double
   computeBalancePenalty(
      const hier::Box& a,
      double imbalance) const
   {
      NULL_USE(a);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

      std::set<BoxInTransit, BoxInTransitMoreLoad> d_set;
      LoadType d_sumload;

   const PartitioningParams *d_pparams;

   BalanceBoxBreaker d_bbb;

   bool d_allow_box_breaking;


   //@{
   //! @name Debugging and diagnostic data.

   bool d_print_steps;
   bool d_print_pop_steps;
   bool d_print_swap_steps;
   bool d_print_break_steps;

   boost::shared_ptr<tbox::Timer> t_adjust_load;
   boost::shared_ptr<tbox::Timer> t_adjust_load_by_popping;
   boost::shared_ptr<tbox::Timer> t_adjust_load_by_swapping;
   boost::shared_ptr<tbox::Timer> t_shift_loads_by_breaking;
   boost::shared_ptr<tbox::Timer> t_find_swap_pair;

   //@}
};



}
}

#endif
