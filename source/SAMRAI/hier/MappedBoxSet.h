/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSet
#define included_hier_MappedBoxSet

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <set>

namespace SAMRAI {
namespace hier {

class BoxList;

/*!
 * @brief A wrapper around std::set<Box>.
 *
 * This is little more than a std::set<Box>, a sorted container
 * of MappedBoxes.  It adds a few additional "features" such as:
 *
 * - Database reading/writing
 * - printing
 */
class MappedBoxSet
{

public:

   /*!
    * @brief Default constructor creates an empty container.
    */
   MappedBoxSet();

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   MappedBoxSet(
      const MappedBoxSet& other);

   //! @brief Destructor.
   virtual ~MappedBoxSet();


   //@{

   //! @name Types defined by std::set.

   typedef std::set<Box, Box::id_less>::iterator iterator;
   typedef std::set<Box, Box::id_less>::const_iterator const_iterator;
   typedef std::set<Box, Box::id_less>::reverse_iterator reverse_iterator;
   typedef std::set<Box, Box::id_less>::const_reverse_iterator const_reverse_iterator;
   typedef std::set<Box, Box::id_less>::key_type key_type;
   typedef std::set<Box, Box::id_less>::value_type value_type;
   typedef std::set<Box, Box::id_less>::size_type size_type;
   typedef std::set<Box, Box::id_less>::reference reference;
   typedef std::set<Box, Box::id_less>::const_reference const_reference;

   //@}


   //@{

   //! @name Set-like interfaces: see STL set documentation.

   /*
    * This is just a subset of the set interface.  Add more as needed.
    * These methods just pass the call off to d_set.
    */

   iterator begin() {
      return d_set.begin();
   }

   iterator end() {
      return d_set.end();
   }

   const_iterator begin() const {
      return d_set.begin();
   }

   const_iterator end() const {
      return d_set.end();
   }

   reverse_iterator rbegin() {
      return d_set.rbegin();
   }

   reverse_iterator rend() {
      return d_set.rend();
   }

   const_reverse_iterator rbegin() const {
      return d_set.rbegin();
   }

   const_reverse_iterator rend() const {
      return d_set.rend();
   }


   iterator insert(
      iterator i,
      const value_type& v) {
      return d_set.insert(i, v);
   }

   bool insert(
      const value_type& v) {
      return d_set.insert(v).second;
   }

   template<class InputIterator>
   void insert(
      InputIterator i,
      InputIterator j) {
      d_set.insert(i, j);
   }

   void erase(
      iterator i) {
      d_set.erase(i);
   }

   size_type erase(
      const key_type& k) {
      return d_set.erase(k);
   }

   void erase(
      iterator first,
      iterator last) {
      d_set.erase(first, last);
   }


   size_t size() const {
      return d_set.size();
   }

   bool empty() const {
      return d_set.empty();
   }

   void clear() {
      d_set.clear();
   }


   iterator find(
      const key_type& k) {
      return d_set.find(k);
   }

   iterator lower_bound(
      const key_type& k) {
      return d_set.lower_bound(k);
   }

   iterator upper_bound(
      const key_type& k) {
      return d_set.upper_bound(k);
   }

   const_iterator find(
      const key_type& k) const {
      return d_set.find(k);
   }


   MappedBoxSet&
   operator = (
      const MappedBoxSet& rhs);

   bool
   operator == (
      const MappedBoxSet& rhs) const;

   bool
   operator != (
      const MappedBoxSet& rhs) const;

   void
   swap(
      MappedBoxSet& other);

   static void
   swap(
      MappedBoxSet& a,
      MappedBoxSet& b);

   //@}


   /*!
    * @brief Whether the subsets of MappedBoxes owned by a given
    * process are the same between this and another MappedBoxSet.
    */
   bool
   isLocallyEqual(
      const MappedBoxSet& other,
      int rank) const;

   /*!
    * @brief Returns the BoxList containing the Boxes from this MappedBoxSet
    * in the requested block.
    */
   tbox::Pointer<BoxList>
   getSingleBlockBoxList(const tbox::Dimension& dim,
      const BlockId& which_block) const;

   /*!
    * @brief Refine the boxes in a MappedBoxSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.  However, if
    * the input and output sets are the same object, replace the
    * entire content with the changed objects.
    *
    * @param[out] output_mapped_boxes
    *
    * @param[in] ratio Ratio in the refinement operation.
    */
   void
   refine(MappedBoxSet& output_mapped_boxes,
      const IntVector& ratio) const;

   /*!
    * @brief Coarsen the boxes in a set<Box>.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.  However, if
    * the input and output sets are the same, replace the entire
    * content with the changed objects.
    *
    * @param[out] output_mapped_boxes
    *
    * @param[in] ratio Ratio in the coarsen operation.
    */
   void
   coarsen(MappedBoxSet& output_mapped_boxes,
      const IntVector& ratio) const;

   /*!
    * @brief Grow the boxes in a MappedBoxSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.  However, if
    * the input and output sets are the same, replace the entire
    * content with the changed objects.
    *
    * @param[out] output_mapped_boxes
    *
    * @param[in] growth Grow boxes by this amount.
    */
   void
   grow(MappedBoxSet& output_mapped_boxes,
      const IntVector& growth) const;

   /*!
    * @brief Unshift periodic image MappedBoxes from a MappedBoxSet.
    *
    * Change periodic image MappedBoxes to their unshifted position.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_mapped_boxes
    *
    * @param[in] input_mapped_boxes
    *
    * @param[in] refinement_ratio Refinement ratio where the boxes
    * live.
    */
   void
   unshiftPeriodicImageMappedBoxes(
      MappedBoxSet& output_mapped_boxes,
      const IntVector& refinement_ratio) const;

   /*!
    * @brief Remove periodic image MappedBoxes from a MappedBoxSet.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_mapped_boxes
    */
   void
   removePeriodicImageMappedBoxes(
      MappedBoxSet& output_mapped_boxes) const;

   /*!
    * @brief Remove from a BoxList its intersections with a MappedBoxSet.
    *
    *
    * @param[in, out] boxes
    */
   void
   removeBoxListIntersections(
      BoxList& boxes) const;

   /*!
    * @brief Insert Box owners into a single set container.
    *
    * @param[out] owners
    */
   void
   getOwners(
      std::set<int>& owners) const;

   /*!
    * @brief Split a MappedBoxSet into two vector<Box>
    * objects, one containing real MappedBoxes and one containing their
    * periodic images.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] real_mapped_box_vector
    *
    * @param[out] periodic_image_mapped_box_vector
    */
   void
   separatePeriodicImages(
      std::vector<Box>& real_mapped_box_vector,
      std::vector<Box>& periodic_image_mapped_box_vector) const;


   //@{
   /*!
    * @name IO support.
    */

   /*!
    * @brief Write the MappedBoxSet to a database.
    */
   void
   putToDatabase(
      tbox::Database& database) const;

   /*!
    * @brief Read the MappedBoxSet from a database.
    */
   void
   getFromDatabase(
      tbox::Database& database);


   /*!
    * @brief Intermediary between MappedBoxSet and output streams,
    * adding ability to control the output.  See
    * MappedBoxSet::format().
    */
   class Outputter {

      friend std::ostream& operator << ( std::ostream& s, const Outputter& f);

   private:

      friend class MappedBoxSet;

      /*!
       * @brief Construct the Outputter with a MappedBoxSet and the
       * parameters needed to output the MappedBoxSet to a stream.
       */
      Outputter( const MappedBoxSet &mapped_box_set,
                 const std::string& border,
                 int detail_depth = 0);

      void operator=( const Outputter &rhs ); // Unimplemented private.

      const MappedBoxSet &d_set;

      const std::string d_border;

      const int d_detail_depth;
   };

   /*!
    * @brief Return a object to that can format the MappedBoxSet for
    * inserting into output streams.
    *
    * Usage example (printing with a tab indentation):
    * @verbatim
    *    cout << "my mapped_boxes:\n" << mapped_boxes.format("\t") << endl;
    * @endverbatim
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   Outputter format( const std::string& border=std::string(),
                     int detail_depth = 0 ) const;

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
      int detail_depth) const;

   //@}

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_MAPPED_BOX_SET_VERSION;

   std::set<Box, Box::id_less> d_set;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxSet.I"
#endif

#endif  // included_hier_MappedBoxSet
