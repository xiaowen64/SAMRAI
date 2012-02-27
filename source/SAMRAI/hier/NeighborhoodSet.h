/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#ifndef included_hier_NeighborhoodSet
#define included_hier_NeighborhoodSet

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/tbox/Database.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <map>
#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A wrapper around std::map<BoxId,BoxContainer>.
 *
 * A neighborhood is defined as a Box and its "neighbors" which
 * are related to it via Connector relationships.  For example,
 * a NeighborhoodSet is used to group the set of Box neighbors
 * for each BoxId (corresponding to a Box).
 *
 * This is little more than a dumb container of neighbor data
 * with a few additional "features" such as:
 * - Database reading/writing
 * - printing
 *
 * The NeighborhoodSet @c m maps a BoxId @c i to a set of Boxes @c
 * m[i], usually its "neighbors".  (The neighbor are stored in a
 * BoxContainer.)
 */

class NeighborhoodSet
{

public:
   //! @brief NeighborSet is a clarifying typedef.
   typedef BoxContainer NeighborSet;

   //! @brief Default constructor creates an empty container.
   NeighborhoodSet();

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   NeighborhoodSet(
      const NeighborhoodSet& other);

   //! @brief Destructor.
   ~NeighborhoodSet();

   //@{

   //! @name Types defined by std::map.

   typedef std::map<BoxId, NeighborSet>::iterator iterator;
   typedef std::map<BoxId, NeighborSet>::const_iterator const_iterator;
   typedef std::map<BoxId, NeighborSet>::reverse_iterator reverse_iterator;
   typedef std::map<BoxId, NeighborSet>::key_type key_type;
   typedef std::map<BoxId, NeighborSet>::value_type value_type;
   typedef std::map<BoxId, NeighborSet>::size_type size_type;

   //@}

   //@{

   //! @name Map-like interfaces: see STL map documentation.

   /*
    * This is just a subset of the map interface.  Add more as needed.
    * These methods just pass the call off to d_map.
    */

   iterator
   begin();

   iterator
   end();

   const_iterator
   begin() const;

   const_iterator
   end() const;

   reverse_iterator
   rbegin();

   reverse_iterator
   rend();

   iterator
   insert(
      iterator i,
      const value_type& v);

   void
   erase(
      iterator i);

   size_type
   erase(
      const key_type& k);

   void
   erase(
      iterator first,
      iterator last);

   size_t
   size() const;

   bool
   empty() const;

   void
   clear();

   NeighborSet&
   operator [] (
      const key_type& k);

   iterator
   find(
      const key_type& k);

   const_iterator
   find(
      const key_type& k) const;

   NeighborhoodSet&
   operator = (
      const NeighborhoodSet& rhs);

   bool
   operator == (
      const NeighborhoodSet& rhs) const;

   bool
   operator != (
      const NeighborhoodSet& rhs) const;

   void
   swap(
      NeighborhoodSet& other);

   //@}

   //@{

   //! @name Structure to represent a range of items in the container.

   /*!
    * @brief Shorthand for an range of neighborhoods in the NeighborhoodSet.
    *
    * A Range has two members, first and second.  it denotes the
    * subset [first,second), which includes first and everything up to
    * but not including second.  An empty range is designated by
    * first == second.
    */
   typedef std::pair<iterator, iterator> Range;

   /*!
    * @brief Shorthand for a range of neighborhoods in a const NeighborhoodSet.
    *
    * A Range has two members, first and second.  it denotes the
    * subset [first,second), which includes first and everything up to
    * but not including second.  An empty range is designated by
    * first == second.
    */
   typedef std::pair<const_iterator, const_iterator> ConstRange;

   //@}

   /*!
    * @brief Find the range of neighborhoods owned by a given rank.
    *
    * The neighborhoods found are in [first,second).  If first==second,
    * no neighborhoods are found for the given rank.
    *
    * @return the range requested.
    */
   Range
   findRanksRange(
      int rank);

   /*!
    * @brief Find the range of neighborhoods owned by a given rank.
    *
    * The neighborhoods found are in [first,second).  If first==second,
    * no neighborhoods are found for the given rank.
    *
    * @return the range requested.
    */
   ConstRange
   findRanksRange(
      int rank) const;

   /*!
    * @brief Refine the boxes of the neighbors in a NeighborhoodSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_neighborhood_set
    *
    * @param[in] ratio Ratio in the refine operation.
    */
   void
   refineNeighbors(
      NeighborhoodSet& output_neighborhood_set,
      const IntVector& ratio) const;

   /*!
    * @brief Coarsen the boxes of the neighbors in a NeighborhoodSet
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_neighborhood_set
    *
    * @param[in] ratio Ratio in the coarsen operation.
    */
   void
   coarsenNeighbors(
      NeighborhoodSet& output_neighborhood_set,
      const IntVector& ratio) const;

   /*!
    * @brief Grow the boxes of the neighbors in a NeighborhoodSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_neighborhood_set
    *
    * @param[in] growth Growth amount.
    */
   void
   growNeighbors(
      NeighborhoodSet& output_neighborhood_set,
      const IntVector& growth) const;

   /*!
    * @brief Remove periodic neighbors from a NeighborhoodSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_neighborhood_set
    */
   void
   removePeriodicNeighbors();

   /*!
    * @brief Insert all neighbors from a NeighborhoodSet into a
    * single NeighborSet.
    *
    * @param[out] all_neighbors
    */
   void
   getNeighbors(
      NeighborSet& all_neighbors) const;

   /*!
    * @brief Insert all neighbors from a NeighborhoodSet with a
    * given BlockId into a single BoxContainer.
    *
    * @param[out] all_neighbors
    */
   void
   getNeighbors(
      BoxContainer& all_neighbors,
      const BlockId& block_id) const;

   /*!
    * @brief Insert all neighbors from a NeighborhoodSet into
    * multiple BoxContainers differentiated by the BlockId of the
    * neighbors.
    *
    * @param[out] all_neighbors
    */
   void
   getNeighbors(
      std::map<BlockId, BoxContainer>& all_neighbors) const;

   /*!
    * @brief Insert all owners of neighbors from a NeighborhoodSet
    * into a single set container.
    *
    * @param[out] owners
    */
   void
   getOwners(
      std::set<int>& owners) const;

   //@{
   /*!
    * @name IO support.
    */

   /*!
    * @brief Write the NeighborhoodSet to a database.
    */
   void
   putUnregisteredToDatabase(
      const boost::shared_ptr<tbox::Database>& database) const;

   /*!
    * @brief Read the NeighborhoodSet from a database.
    */
   void
   getFromDatabase(
      tbox::Database& database);

   //@}

   /*!
    * @brief Intermediary between NeighborhoodSet and output streams,
    * adding ability to control the output.  See
    * NeighborhoodSet::format().
    */
   class Outputter
   {

      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);

private:
      friend class NeighborhoodSet;

      /*!
       * @brief Construct the Outputter with a NeighborhoodSet and the
       * parameters needed to output the NeighborhoodSet to a stream.
       */
      Outputter(
         const NeighborhoodSet& neighborhood_set,
         const std::string& border,
         int detail_depth = 0);

      void
      operator = (
         const Outputter& rhs);               // Unimplemented private.

      const NeighborhoodSet& d_set;

      const std::string d_border;

      const int d_detail_depth;
   };

   /*!
    * @brief Return a object to that can format the NeighborhoodSet for
    * inserting into output streams.
    *
    * Usage example (printing with a tab indentation):
    * @verbatim
    *    cout << "my neighborhoods:\n" << neighborhoods.format("\t") << endl;
    * @endverbatim
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 0) const;

   /*!
    * @brief Print the contents of the object recursively.
    *
    * @param[in] output_stream
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   void
   recursivePrint(
      std::ostream& output_stream,
      const std::string& left_border,
      int detail_depth = 0) const;

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_EDGE_SET_VERSION;

   /*!
    * @brief NeighborhoodSet is just a wrapper around an STL map.
    */
   std::map<BoxId, NeighborSet> d_map;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/NeighborhoodSet.I"
#endif

#endif // included_hier_NeighborhoodSet
