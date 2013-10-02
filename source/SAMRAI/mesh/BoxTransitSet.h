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
#include "SAMRAI/mesh/TreeLoadBalancer.h"

#include <set>

namespace SAMRAI {
namespace mesh {


/*
 * @brief Implementation of TreeLoadBalancer::TransitSet, representing
 * the load with a set of boxes.
 */

class BoxTransitSet {

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
      BoxInTransit(const tbox::Dimension& dim);

      /*!
       * @brief Construct a new BoxInTransit from an originating box.
       *
       * @param[in] origin
       */
      BoxInTransit(const hier::Box& origin);

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
         hier::LocalId local_id);

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
      putToMessageStream(
         tbox::MessageStream &msg) const;

      /*!
       * @brief Set attributes according to data in a MessageStream.
       *
       * This is the opposite of putToMessageStream().
       */
      void
      getFromMessageStream(
         tbox::MessageStream &msg);

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
      const BoxInTransit& r);


   /*!
    * @brief Comparison functor for sorting BoxInTransit from more to
    * less loads.
    */
   struct BoxInTransitMoreLoad {
      /*
       * @brief Compares two BoxInTransit for sorting them from more load
       * to less load.
       */
      bool
      operator () (
         const BoxInTransit& a,
         const BoxInTransit& b) const
      {
         if (a.getBox().size() != b.getBox().size()) {
            return a.d_boxload > b.d_boxload;
         }
         return a.d_box.getBoxId() < b.d_box.getBoxId();
      }
   };




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
      BoxTransitSet() : d_set(), d_sumload(0) {}
      template<class InputIterator>
      BoxTransitSet( InputIterator first, InputIterator last ) :
         d_set(first,last), d_sumload(0) {
         for ( const_iterator bi=d_set.begin(); bi!=d_set.end(); ++bi )
         { d_sumload += bi->d_boxload; };
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
   private:
      std::set<BoxInTransit, BoxInTransitMoreLoad> d_set;
      LoadType d_sumload;
   };



}
}

#endif
