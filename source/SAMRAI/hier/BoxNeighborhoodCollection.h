/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A class describing the adjacency of Boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxNeighborhoodCollection
#define included_hier_BoxNeighborhoodCollection

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/BoxSet.h"
#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <map>
#include <set>
#include <vector>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Given a box called the root, the boxes which are adjacent to the root
 * are its neighbors and are said to form the neighborhood of the root.  The
 * neighborhood of a root Box does not include the root Box itself.  This class
 * describes the neighborhoods of a collection of root Boxes.  Each root in the
 * collection has a neighborhood of adjacent Boxes.
 */
class BoxNeighborhoodCollection
{
   friend class Iterator;
   friend class ConstIterator;

   private:
      // Default constructor does not exist.
      BoxNeighborhoodCollection();

      /*!
       * @brief Rank of process on which this object lives.
       */
      int d_rank;

      // Hmmmm.
      struct box_ptr_less {
         bool operator () (const Box* box0, const Box* box1) const
         {
            return box0->getId() < box1->getId();
         }
      };

      struct box_less {
         bool operator () (const Box& box0, const Box& box1) const
         {
            return box0.getId() < box1.getId();
         }
      };

      // Typedefs.

      typedef std::set<const Box*, box_ptr_less> Roots;

      typedef Roots::iterator RootsItr;

      typedef std::set<const Box*, box_ptr_less> Neighbors;

      typedef Neighbors::iterator NeighborsItr;

      typedef Neighbors::const_iterator NeighborsConstItr;

      typedef std::map<Box, Neighbors, box_less> AdjList;

      typedef AdjList::iterator AdjListItr;

      typedef AdjList::const_iterator AdjListConstItr;

      AdjList d_adj_list;

      Roots d_roots;

   public:
      // Constructors.

      /*!
       * @brief Constructs an empty object. There are not yet any root Boxes
       * whose neighborhoods are represented by this object.
       *
       * @param rank The rank of the processor on which this object lives.
       */
      BoxNeighborhoodCollection(
         int rank);

      /*!
       * @brief Constructs a collection of empty neighborhoods for each root
       * box in roots.
       *
       * @param rank The rank of the processor on which this object lives.
       *
       * @param roots Root boxes whose neighborhoods will be represented by
       * this object.
       */
      BoxNeighborhoodCollection(
         int rank,
         const BoxContainer& roots);

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      BoxNeighborhoodCollection(
         const BoxNeighborhoodCollection& other);

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      BoxNeighborhoodCollection& operator = (
         const BoxNeighborhoodCollection& rhs);


      // Destructor

      /*!
       * @brief Destructor
       */
      ~BoxNeighborhoodCollection();


      // Operators

      /*!
       * @brief Determine if two collections are equivalent.
       *
       * @param rhs
       */
      bool operator == (
         const BoxNeighborhoodCollection& rhs) const;


      // Iteration

      /*!
       * @brief An iterator over the roots of the neighborhoods in a
       * BoxNeighborhoodCollection.  The interface does not allow modification
       * of the neighborhood roots.
       */
      class Iterator
      {
         friend class BoxNeighborhoodCollection;

         public:
            // Constructors.

            /*!
             * @brief Constructs an iterator over the neighborhood roots in
             * the supplied collection.
             *
             * @param nbrhds
             */
            Iterator(
               BoxNeighborhoodCollection& nbrhds,
               bool from_start = true);

            /*!
             * @brief Copy constructor.
             *
             * @param other
             */
            Iterator(
               const Iterator& other);

            /*!
             * @brief Assignment operator.
             *
             * @param rhs
             */
            Iterator&
            operator = (
               const Iterator& rhs);


            // Destructor

            /*!
             * @brief Performs necessary deletion.
             */
            ~Iterator();


            // Operators

            /*!
             * @brief Extracts the Box which is the root of the current
             * neighborhood in the iteration.
             */
            const Box&
            operator * ();

            /*!
             * @brief Pre-increment iterator to point to Box which is root of
             * next neighborhood in the collection.
             */
            Iterator&
            operator ++ ();

            /*!
             * @brief Determine if two iterators are equivalent.
             *
             * @param rhs
             */
            bool
            operator == (
               const Iterator& rhs) const;

            /*!
             * @brief Determine if two iterators are not equivalent.
             *
             * @param rhs
             */
            bool
            operator != (
               const Iterator& rhs) const;

         private:
            // Default constructor does not exist.
            Iterator();

            /*!
             * @brief Constructs an iterator pointing to a specific root in
             * nbrhds.  Should only be called by BoxNeighborhoodCollection.
             *
             * @param nbrhds
             *
             * @param itr
             */
            Iterator(
               BoxNeighborhoodCollection& nbrhds,
               AdjListItr itr);

            const BoxNeighborhoodCollection* d_collection;

            AdjListItr d_itr;
      };

      /*!
       * @brief An iterator over the roots of the neighborhoods in a const
       * BoxNeighborhoodCollection.  The interface does not allow modification
       * of the neighborhood roots.
       */
      class ConstIterator
      {
         friend class BoxNeighborhoodCollection;

         public:
            // Constructors.

            /*!
             * @brief Constructs an iterator over the neighborhood roots in
             * the supplied collection.
             *
             * @param nbrhds
             */
            ConstIterator(
               const BoxNeighborhoodCollection& nbrhds,
               bool from_start = true);

            /*!
             * @brief Copy constructor.
             *
             * @param other
             */
            ConstIterator(
               const ConstIterator& other);

            /*!
             * @brief Assignment operator.
             *
             * @param rhs
             */
            ConstIterator&
            operator = (
               const ConstIterator& rhs);


            // Destructor

            /*!
             * @brief Performs necessary deletion.
             */
            ~ConstIterator();


            // Operators

            /*!
             * @brief Extracts the Box which is the root of the current
             * neighborhood in the iteration.
             */
            const Box&
            operator * () const;

            /*!
             * @brief Pre-increment iterator to point to Box which is root of
             * next neighborhood in the collection.
             */
            ConstIterator&
            operator ++ ();

            /*!
             * @brief Determine if two iterators are equivalent.
             *
             * @param rhs
             */
            bool
            operator == (
               const ConstIterator& rhs) const;

            /*!
             * @brief Determine if two iterators are not equivalent.
             *
             * @param rhs
             */
            bool
            operator != (
               const ConstIterator& rhs) const;

         private:
            // Default constructor does not exist.
            ConstIterator();

            /*!
             * @brief Constructs an iterator pointing to a specific root in
             * nbrhds.  Should only be called by BoxNeighborhoodCollection.
             *
             * @param nbrhds
             *
             * @param itr
             */
            ConstIterator(
               const BoxNeighborhoodCollection& nbrhds,
               AdjListConstItr itr);

            const BoxNeighborhoodCollection* d_collection;

            AdjListConstItr d_itr;
      };

      /*!
       * @brief An iterator over the neighbors in the neighborhood of a root in
       * a BoxNeighborhoodCollection.  The interface does not allow
       * modification of the neighbors.
       */
      class NeighborIterator
      {
         friend class BoxNeighborhoodCollection;

         public:
            // Constructors

            /*!
             * @brief Constructs an iterator over the neighbors of the supplied
             * root in the supplied collection of neighborhoods.
             *
             * @param nbrhds
             *
             * @param root
             */
            NeighborIterator(
               BoxNeighborhoodCollection& nbrhds,
               Iterator& root,
               bool from_start = true);

            /*!
             * @brief Copy constructor.
             *
             * @param other
             */
            NeighborIterator(
               const NeighborIterator& other);

            /*!
             * @brief Assignment operator.
             *
             * @param rhs
             */
            NeighborIterator&
            operator = (
               const NeighborIterator& rhs);


            // Destructor

            /*!
             * @brief Performs necessary deletion.
             */
            ~NeighborIterator();


            // Operators

            /*!
             * @brief Extract the Box which is the current neighbor in the
             * iteration of the neighborhood of the root.
             */
            const Box&
            operator * () const;

            /*!
             * @brief Pre-increment iterator to point to the Box which is the
             * next neighbor of the root.
             */
            NeighborIterator&
            operator ++ ();

            /*!
             * @brief Determine if two iterators are equivalent.
             *
             * @param rhs
             */
            bool
            operator == (
               const NeighborIterator& rhs) const;

            /*!
             * @brief Determine if two iterators are not equivalent.
             *
             * @param rhs
             */
            bool
            operator != (
               const NeighborIterator& rhs) const;

         private:
            // Default constructor does not exist.
            NeighborIterator();

            const BoxNeighborhoodCollection* d_collection;

            const Box* d_root;

            NeighborsItr d_itr;
      };

      /*!
       * @brief An iterator over the neighbors in the neighborhood of a root in
       * a const BoxNeighborhoodCollection.  The interface does not allow
       * modification of the neighbors.
       */
      class ConstNeighborIterator
      {
         friend class BoxNeighborhoodCollection;

         public:
            // Constructors

            /*!
             * @brief Constructs an iterator over the neighbors of the supplied
             * root in the supplied collection of neighborhoods.
             *
             * @param nbrhds
             *
             * @param root
             */
            ConstNeighborIterator(
               const BoxNeighborhoodCollection& nbrhds,
               const ConstIterator& root,
               bool from_start = true);

            /*!
             * @brief Copy constructor.
             *
             * @param other
             */
            ConstNeighborIterator(
               const ConstNeighborIterator& other);

            /*!
             * @brief Assignment operator.
             *
             * @param rhs
             */
            ConstNeighborIterator&
            operator = (
               const ConstNeighborIterator& rhs);


            // Destructor

            /*!
             * @brief Performs necessary deletion.
             */
            ~ConstNeighborIterator();


            // Operators

            /*!
             * @brief Extract the Box which is the current neighbor in the
             * iteration of the neighborhood of the root.
             */
            const Box&
            operator * () const;

            /*!
             * @brief Pre-increment iterator to point to the Box which is the
             * next neighbor of the root.
             */
            ConstNeighborIterator&
            operator ++ ();

            /*!
             * @brief Determine if two iterators are equivalent.
             *
             * @param rhs
             */
            bool
            operator == (
               const ConstNeighborIterator& rhs) const;

            /*!
             * @brief Determine if two iterators are not equivalent.
             *
             * @param rhs
             */
            bool
            operator != (
               const ConstNeighborIterator& rhs) const;

         private:
            // Default constructor does not exist.
            ConstNeighborIterator();

            const BoxNeighborhoodCollection* d_collection;

            const Box* d_root;

            NeighborsConstItr d_itr;
      };

      /*!
       * @brief Returns an iterator pointing to the beginning of the collection
       * of neighborhoods.
       */
      Iterator
      begin();

      /*!
       * @brief Returns an iterator pointing to the beginning of the collection
       * of neighborhoods.
       */
      ConstIterator
      begin() const;

      /*!
       * @brief Returns an iterator pointing just past the end of the
       * collection of neighborhoods.
       */
      Iterator
      end();

      /*!
       * @brief Returns an iterator pointing just past the end of the
       * collection of neighborhoods.
       */
      ConstIterator
      end() const;

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of root.
       *
       * @param root
       */
      NeighborIterator
      begin(
         const Box& root);

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of root.
       *
       * @param root
       */
      ConstNeighborIterator
      begin(
         const Box& root) const;

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of the root Box with the supplied BoxId.
       *
       * @param root_id
       */
      NeighborIterator
      begin(
         const BoxId& root_id);

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of the root Box with the supplied BoxId.
       *
       * @param root_id
       */
      ConstNeighborIterator
      begin(
         const BoxId& root_id) const;

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of the root Box pointed to by root_itr.
       *
       * @param root_itr
       */
      NeighborIterator
      begin(
         Iterator& root_itr);

      /*!
       * @brief Returns an iterator pointing to the first neighbor in the
       * neighborhood of the root Box pointed to by root_itr.
       *
       * @param root_itr
       */
      ConstNeighborIterator
      begin(
         const ConstIterator& root_itr) const;

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of root.
       *
       * @param root
       */
      NeighborIterator
      end(
         const Box& root);

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of root.
       *
       * @param root
       */
      ConstNeighborIterator
      end(
         const Box& root) const;

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of the root Box with the supplied BoxId.
       *
       * @param root_id
       */
      NeighborIterator
      end(
         const BoxId& root_id);

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of the root Box with the supplied BoxId.
       *
       * @param root_id
       */
      ConstNeighborIterator
      end(
         const BoxId& root_id) const;

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of the root Box pointed to by root_itr.
       *
       * @param root_itr
       */
      NeighborIterator
      end(
         Iterator& root_itr);

      /*!
       * @brief Returns an iterator pointing just past the last neighbor in the
       * neighborhood of the root Box pointed to by root_itr.
       *
       * @param root_itr
       */
      ConstNeighborIterator
      end(
         const ConstIterator& root_itr) const;


      // Lookup

      /*!
       * @brief Returns an iterator pointing to the supplied neighborhood root.
       * If root is not in the collection this method returns end().
       *
       * @param root
       */
      ConstIterator
      find(
         const Box& root) const;

      /*!
       * @brief Returns an iterator pointing to the neighborhood root with the
       * supplied BoxId.  If no root's BoxId is root_id this method returns
       * end().
       *
       * @param root_id
       */
      ConstIterator
      find(
         const BoxId& root_id) const;

      /*!
       * @brief Returns an iterator pointing to the supplied neighborhood root.
       * If root is not in the collection this method returns end().
       *
       * @param root
       */
      Iterator
      find(
         const Box& root);

      /*!
       * @brief Returns an iterator pointing to the neighborhood root with the
       * supplied BoxId.  If no root's BoxId is root_id this method returns
       * end().
       *
       * @param root_id
       */
      Iterator
      find(
         const BoxId& root_id);


      // Typedefs
      typedef std::pair<Iterator, bool> InsertRetType;


      // State queries

      /*!
       * @brief Returns true if the number of box neighborhoods == 0.
       */
      bool
      empty() const;

      /*!
       * @brief Returns the number of box neighborhoods.
       */
      int
      numBoxNeighborhoods() const;

      /*!
       * @brief Returns true if the neighborhood of root is empty.
       *
       * @param root
       */
      bool
      emptyBoxNeighborhood(
         const Box& root) const;

      /*!
       * @brief Returns true if the neighborhood of the root Box with the
       * supplied BoxId is empty.
       *
       * @param root_id
       */
      bool
      emptyBoxNeighborhood(
         const BoxId& root_id) const;

      /*!
       * @brief Returns true if the neighborhood of the root Box pointed to by
       * root_itr is empty.
       *
       * @param root_itr
       */
      bool
      emptyBoxNeighborhood(
         const ConstIterator& root_itr) const;

      /*!
       * @brief Returns the number of neighbors in the neighborhood of root.
       *
       * @param root
       */
      int
      numNeighbors(
         const Box& root) const;

      /*!
       * @brief Returns the number of neighbors in the neighborhood of the root
       * Box with the supplied BoxId.
       *
       * @param root_id
       */
      int
      numNeighbors(
         const BoxId& root_id) const;

      /*!
       * @brief Returns the number of neighbors in the neighborhood of the root
       * Box pointed to by root_itr.
       *
       * @param root_itr
       */
      int
      numNeighbors(
         const ConstIterator& root_itr) const;

      /*!
       * @brief Returns the number of neighbors in all neighborhoods.
       */
      int
      totalNumNeighbors() const;

      /*!
       * @brief Returns true if any neighbor of root is a periodic Box.
       *
       * @param root
       */
      bool
      hasPeriodicNeighborhood(
         const Box& root) const;

      /*!
       * @brief Returns true if any neighbor of the root Box with the supplied
       * BoxId is a periodic Box.
       *
       * @param root_id
       */
      bool
      hasPeriodicNeighborhood(
         const BoxId& root_id) const;

      /*!
       * @brief Returns true if any neighbor of root Box pointed to by root_itr
       * is a periodic Box.
       *
       * @param root_itr
       */
      bool
      hasPeriodicNeighborhood(
         const ConstIterator& root_itr) const;

      /*!
       * @brief Returns true if root is the root of a neighborhood held by this
       * object.
       *
       * @param root
       */
      bool
      isRoot(
         const Box& root) const;

      /*!
       * @brief Returns true if root is the BoxId of the root of a neighborhood
       * held by this object.
       *
       * @param root_id
       */
      bool
      isRoot(
         const BoxId& root_id) const;

      /*!
       * @brief Returns true if all neighbors of all roots are owned by the
       * processor owning this object.
       */
      bool
      isLocal() const;

      /*!
       * @brief Insert the rank of the processor owning each neighbor in each
       * neighborhood into the supplied set.
       *
       * @param owners
       */
      void
      getOwners(
         std::set<int>& owners) const;


      // Neighborhood editing

      /*!
       * @brief Inserts a new neighbor into the neighborhood of root.  If the
       * neighborhood of root is not yet represented by this object, then the
       * object will be modified to include this new box neighborhood and the
       * method will return true.  Otherwise it will return false.
       *
       * @param root The neighborhood root.
       *
       * @param new_nbr The new neighbor of root.
       *
       * @return A pair the first member of which is an Iterator pointing to
       * root and the second of which is true if root is a new neighborhood
       * root, and false otherwise.
       */
      InsertRetType
      insert(
         const Box& root,
         const Box& new_nbr);

      /*!
       * @brief Inserts a new neighbor into the neighborhood of the root with
       * the supplied BoxId.
       *
       * @note root_id must be the BoxId of a valid root or this function can
       * not work.  We can't create a new root Box only knowing a BoxId.
       * Therefore this version of insert only returns an Iterator pointing to
       * the root.
       *
       * @param root_id The BoxId of the neighborhood root.
       *
       * @param new_nbr The new neighbor of root.
       *
       * @return An Iterator pointing to the root with the supplied BoxId.
       */
      Iterator
      insert(
         const BoxId& root_id,
         const Box& new_nbr);

      /*!
       * @brief Inserts a new neighbor into the neighborhood of the root
       * pointed to by root_itr.
       *
       * @note root_itr must point to a valid root or this function can not
       * work.  Unlike the other versions of insert, this version has no return
       * value.  The root must exist and the Iterator already points to it so
       * returning an Iterator or bool has no value.
       *
       * @param root_itr Iterator pointing to the neighborhood root.
       *
       * @param new_nbr The new neighbor of root.
       */
      void
      insert(
         Iterator& root_itr,
         const Box& new_nbr);

      /*!
       * @brief Inserts new neighbors into the neighborhood of root.  If the
       * neighborhood of root is not yet represented by this object, then the
       * object will be modified to include this new box neighborhood and the
       * method will return true.  Otherwise it will return false.
       *
       * @param root The neighborhood root.
       *
       * @param new_nbrs The new neighbors of root.
       *
       * @return A pair the first member of which is an Iterator pointing to
       * root and the second of which is true if root is a new neighborhood
       * root, and false otherwise.
       */
      InsertRetType
      insert(
         const Box& root,
         const BoxContainer& new_nbrs);

      /*!
       * @brief Inserts new neighbors into the neighborhood of the root with
       * the supplied BoxId.
       *
       * @note root_id must be the BoxId of a valid root or this function can
       * not work.  We can't create a new root Box only knowing a BoxId.
       * Therefore this version of insert only returns an Iterator pointing to
       * the root.
       *
       * @param root_id The BoxId of the neighborhood root.
       *
       * @param new_nbrs The new neighbors of root.
       *
       * @return An Iterator pointing to the root with the supplied BoxId.
       */
      Iterator
      insert(
         const BoxId& root_id,
         const BoxContainer& new_nbrs);

      /*!
       * @brief Inserts new neighbors into the neighborhood of the root pointed
       * to by root_itr.
       *
       * @note root_itr must point to a valid root or this function can not
       * work.  Unlike the other versions of insert, this version has no return
       * value.  The root must exist and the Iterator already points to it so
       * returning an Iterator or bool has no value.
       *
       * @param root_itr Iterator pointing to the neighborhood root.
       *
       * @param new_nbrs The new neighbors of root.
       */
      void
      insert(
         Iterator& root_itr,
         const BoxContainer& new_nbrs);

      /*!
       * @brief Erases a neighbor from the neighborhood of root.
       *
       * @param root The neighborhood root.
       *
       * @param nbr The neighbor of root to be erased.
       */
      void
      erase(
         const Box& root,
         const Box& nbr);

      /*!
       * @brief Erases a neighbor from the neighborhood of the root with the
       * supplied BoxId.
       *
       * @param root_id The BoxId of the neighborhood root.
       *
       * @param nbr The neighbor of the root to be erased.
       */
      void
      erase(
         const BoxId& root_id,
         const Box& nbr);

      /*!
       * @brief Erases a neighbor from the neighborhood of the root pointed
       * to by root_itr.
       *
       * @param root_itr An iterator pointing to the neighborhood root.
       *
       * @param nbr The neighbor of root to be erased.
       */
      void
      erase(
         Iterator& root_itr,
         const Box& nbr);

      /*!
       * @brief Erases neighbors from the neighborhood of root.
       *
       * @param root The neighborhood root.
       *
       * @param nbrs The neighbors of root to be erased.
       */
      void
      erase(
         const Box& root,
         const BoxContainer& nbrs);

      /*!
       * @brief Erases neighbors from the neighborhood of the root with the
       * supplied BoxId.
       *
       * @param root_id The BoxId of the neighborhood root.
       *
       * @param nbrs The neighbors of root to be erased.
       */
      void
      erase(
         const BoxId& root_id,
         const BoxContainer& nbrs);

      /*!
       * @brief Erases neighbors from the neighborhood of the root pointed to
       * by root_itr.
       *
       * @param root_itr An iterator pointing to the neighborhood root.
       *
       * @param nbrs The neighbors of root to be erased.
       */
      void
      erase(
         Iterator& root_itr,
         const BoxContainer& nbrs);

      /*!
       * @brief Inserts a root with an empty neighborhood.  If the root does
       * not exist in this object then this function returns true.  Otherwise
       * it returns false and takes no action.
       *
       * @param new_root The new root of an empty neighborhood.
       *
       * @return A pair the first member of which is an Iterator pointing to
       * root and the second of which is true if root is a new neighborhood
       * root, and false otherwise.
       */
      InsertRetType
      insert(
         const Box& new_root);

      /*!
       * @brief Erases the neighbors of root.  If erase_root is true then the
       * root will be erased as well.
       *
       * @param root The neighborhood root whose neighbors are to be erased.
       *
       * @param erase_root
       */
      void
      erase(
         const Box& root,
         bool erase_root);

      /*!
       * @brief Erases the neighbors of the root with the supplied BoxId.  If
       * erase_root is true then the root will be erased as well.
       *
       * @param root_id The BoxId of the neighborhood root whose neighbors are
       * to be erased.
       *
       * @param erase_root
       */
      void
      erase(
         const BoxId& root_id,
         bool erase_root);

      /*!
       * @brief Erases the neighbors of the root pointed to by root_itr.  If
       * erase_root is true then the root will be erased as well and root will
       * point to end().
       *
       * @param root_itr Iterator pointing to the neighborhood root whose
       * neighbors are to be erased.
       *
       * @param erase_root
       */
      void
      erase(
         Iterator& root_itr,
         bool erase_root);

      /*!
       * @brief Erases the neighbors of the roots pointed to in the range
       * [first_root_itr, last_root_itr).  If erase_roots is true then the
       * roots will be erased as well and first_root_itr will point to end().
       *
       * @param first_root_itr Iterator pointing to the first root whose
       * neighbors are to be erased.
       *
       * @param last_root_itr Iterator pointing one past the last root whose
       * neighbors are to be erased.
       *
       * @param erase_roots
       */
      void
      erase(
         Iterator& first_root_itr,
         Iterator& last_root_itr,
         bool erase_roots);

      /*!
       * @brief For all roots not owned by the process owning this object this
       * method erases the root and its neighbors.
       */
      void
      localize();

      /*!
       * @brief Erases all roots having no neighbors.
       */
      void
      eraseEmptyNeighborhoods();

      /*!
       * @brief Erases all contents so empty() == true.
       */
      void
      clear();


      // Coarsen, refine, grow.

      /*!
       * @brief Produces another BoxNeighborhoodCollection by coarsening the
       * neighbors of each root by the given ratio.
       *
       * @param result
       *
       * @param ratio
       */
      void
      coarsenNeighbors(
         BoxNeighborhoodCollection& result,
         const IntVector& ratio) const;

      /*!
       * @brief Produces another BoxNeighborhoodCollection by coarsening the
       * neighbors of each root by the given ratio.
       *
       * @param result
       *
       * @param ratio
       */
      void
      refineNeighbors(
         BoxNeighborhoodCollection& result,
         const IntVector& ratio) const;

      /*!
       * @brief Produces another BoxNeighborhoodCollection by growing the
       * neighbors of each root by the given amount.
       *
       * @param result
       *
       * @param growth
       */
      void
      growNeighbors(
         BoxNeighborhoodCollection& result,
         const IntVector& growth) const;


      // Neighborhood member extraction
      // Currently, the way some algorithms are implemented these are needed
      // but we may find that it is not necessary.

      /*!
       * @brief Fill the supplied BoxSet with the neighbors from all the
       * neighborhoods in this object.  The container has no notion of the
       * neighborhoods to which its contents belong.
       *
       * @param neighbors
       */
      void
      getNeighbors(
         BoxSet& neighbors) const;

      /*!
       * @brief Fill the supplied BoxList with the neighbors from all the
       * neighborhoods in this object.  The container has no notion of the
       * neighborhoods to which its contents belong.
       *
       * @param neighbors
       */
//      void
//      getNeighbors(
//         BoxList& neighbors) const;

      /*!
       * @brief Fill the supplied BoxList with the neighbors having the
       * specified block id from all the neighborhoods in this object.  The
       * container has no notion of the neighborhoods to which its contents
       * belong.
       *
       * @param neighbors
       *
       * @param block_id
       */
      void
      getNeighbors(
         BoxList& neighbors,
         const BlockId& block_id) const;

      /*!
       * @brief Fill the supplied map with the neighbors from all the
       * neighborhoods in this object by block id.  The container has no notion
       * of the neighborhoods to which its contents belong.
       *
       * @param neighbors
       */
      void
      getNeighbors(
         std::map<BlockId, BoxList>& neighbors) const;

      /*!
       * @brief Place any periodic neighbors from each neighborhood into the
       * supplied BoxSet.
       *
       * @param result
       */
      void
      getPeriodicNeighbors(
         BoxSet& result) const;


      // Communication packing/unpacking

      /*!
       * @brief Load an integer communication buffer with the data from this
       * object.
       *
       * @param send_mesg The integer communication buffer
       *
       * @param dim
       *
       * @param offset Starting location in send_mesg
       *
       * @param buff_init Initializer for newly allocated buffer data.
       */
      void
      putToIntBuffer(
         std::vector<int>& send_mesg,
         const tbox::Dimension& dim,
         int offset,
         int buff_init) const;

      /*!
       * @brief Populate object based on information contained in an integer
       * communication buffer.
       *
       * @param recv_mesg The integer communication buffer
       *
       * @param dim
       */
      void
      getFromIntBuffer(
         const std::vector<int>& recv_mesg,
         const std::vector<int>& proc_offset,
         const tbox::Dimension& dim,
         int num_proc);


      // IO
      // These are defined in NeighborhoodSet but are only called from the
      // mblktree test.

      /*!
       * @brief Writes the neighborhood information to the supplied database.
       *
       * @param database
       */
      void
      putToDatabase(
         tbox::Database& database) const;

      /*!
       * @brief Constructs the neighborhoods from the supplied database.
       *
       * @param database
       */
      void
      getFromDatabase(
         const tbox::Database& database);
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxNeighborhoodCollection.I"
#endif

#endif // included_hier_BoxNeighborhoodCollection
