/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence classes. 
 *
 ************************************************************************/

#ifndef included_xfer_RefineClasses
#define included_xfer_RefineClasses

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <iostream>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Maintain a collection of refine items and organize them
 * into equivalence classes.
 *
 * RefineClasses is used by the RefineSchedule and RefineAlgorithm
 * classes to manage refinement data items that describe communication 
 * of patch data on an AMR hierarchy.  Specifically, this class organizes
 * these items into equivalence clases, so that items are grouped
 * together if they are considered equivalent.
 *
 * Two items are equivalent if all of the following are true:
 * <ul>
 *   <li> The referenced patch data type for the source component
 *        is the same for both items.
 *   <li> The referenced patch data type for the destination component
 *        is the same for both items.
 *   <li> The ghost cell width of the source patch data is the same
 *        for both items.
 *   <li> The ghost cell width of the destination patch data is the same
 *        for both items.
 *   <li> They use the same VariableFillPattern.
 * </ul>
 */

class RefineClasses:public tbox::DescribedClass
{
public:
   /*!
    * @brief Data structure used to describe a refinement operation
    * between patch data components on an AMR hierarchy.
    */
   struct Data {
      /*!
       * @brief Destination patch data component
       */
      int d_dst;

      /*!
       * @brief Source patch data component
       */
      int d_src;

      /*!
       * @brief Patch data component for source data at the old time in
       * a time interpolation operation.
       */
      int d_src_told;

      /*!
       * @brief Patch data component for source data at the new time in
       * a time interpolation operation.
       */
      int d_src_tnew;

      /*!
       * @breif Scratch patch data component
       */
      int d_scratch;

      /*!
       * @brief Boolean flag that is set to true when it is desired that fine
       * data values have priority over coarse data values when data exists
       * at the same location in the mesh on levels with different resolutions
       * (e.g., nodes that are coincident on consecutive mesh levels).
       */
      bool d_fine_bdry_reps_var;

      /*!
       * @brief Boolean flag telling if this item uses time interpolation.
       */
      bool d_time_interpolate;

      /*!
       * @brief Refinement operator
       */
      tbox::Pointer<hier::RefineOperator> d_oprefine;

      /*!
       * @brief Time interpolation operator
       */
      tbox::Pointer<hier::TimeInterpolateOperator> d_optime;

      /*!
       * @brief Identifier of equivalence class where this item belongs.  All
       * items of the same equivalence class will have the same value.
       */ 
      int d_class_id;

      /*!
       * @brief An array index telling where this item sits in an array of
       * refine items.
       */
      int d_tag;

      /*!
       * @brief VariableFillPattern that can restrict the stencil of the data
       * filled by the RefineSchedule. 
       */ 
      tbox::Pointer<VariableFillPattern> d_var_fill_pattern;
   };

   /*!
    * @brief The constructor creates an empty array of refine classes.
    */
   RefineClasses();

   /*!
    * @brief The virtual destructor destroys the refinement data items owned
    * by this object.
    */
   virtual ~RefineClasses();

   /*!
    * @brief Return number of equivalence classes maintained by this object.
    */
   int
   getNumberOfEquivalenceClasses() const;

   /*!
    * @brief Return total number of refine items that have been registered and
    * stored in the RefineClasses object
    */
   int
   getNumberOfRefineItems() const;

   /*!
    * @brief Get a representative item for a given equivalence class.
    *
    * When assertion checking is active, the id will be checked for validity.
    *
    * @return Given an id indicating a specific equivalence class, one item
    * from that class is returned, to represent the characteristics of the
    * equivalence class.
    *
    * @param[in] equiv_class_id
    */
   const RefineClasses::Data&
   getClassRepresentative(
      int equiv_class_id) const;

   /*!
    * @brief Get a refine item from the array of all refine items held by
    * this object.
    *
    * The internal storage of the refine items held by this class is not
    * controlled by the user, so this method is intended for use when looping
    * over all of the items, from 0 to getNumberOfRefineItems()-1, or when
    * looping over the integers in the List obtained from getIterator().
    *
    * @return A refine item identified by an integer id.
    *
    * @param[in] refine_item_array_id
    */
   RefineClasses::Data&
   getRefineItem(
      const int refine_item_array_id);

   /*!
    * @brief Return an iterator for the list of array ids corresponding to the
    * equivalence class with the given integer identifier.
    *
    * The number of quivalence classes can be determined via the
    * getNumberOfEquivalenceClasses() member function.  Valid integer
    * arguments are from 0 to getNumberOfEquivalenceClasses()-1.  When
    * assertion checking is active, the id will be checked for validity.
    *
    * @note The list should not be modified through this iterator.
    *
    * @return The iterator iterates over a list of integers which are array
    * ids that can be passed into getRefineItem().  The array ids in a
    * single list all correspond to refine items in a single equivalence
    * class.
    *
    * @param[in] equiv_class_id
    */
   tbox::List<int>::Iterator
   getIterator(
      int equiv_class_id);

   /*!
    * @brief Give a data item to the RefineClasses object, which will store
    * it with the proper equivalence class.
    *
    * If the item belongs in an existing equivalence class, it will be added
    * there, otherwise a new equivalence class will be created for this item.
    * The internal data of the item will be changed so that it stores an
    * integer identifier of its equivalence class.
    *
    * An error will occur with a descriptive message if the data item is
    * not valid.  See checkRefineItem() for explanation of validity.
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will be
    * used.
    *
    * @param[in,out] data
    * @param[in] descriptor
    */
   void
   insertEquivalenceClassItem(
      RefineClasses::Data& data_item,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL));

   /*!
    * @brief Check refine data item for validity.
    *
    * A refine data item is invalid if any of its patch data components are
    * negative, or if its scratch data does not have sufficient ghost width
    * for the stencil of the refine operator or the fill pattern, or if the
    * data types of the source and destination data are not compatible
    * to be copied from one to another.
    *
    * An error will occur with a descriptive message if the item is invalid.
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return True if the item is valid.
    *
    * @param[in] data_item
    * @param[in] descriptor
    */
   bool
   checkRefineItem(
      const RefineClasses::Data& data_item,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL)) const;

   /*!
    * @brief Compare RefineClasses object with another RefineClasses object.
    *
    * This method checks if the equivalence classes held by the two objects
    * match with regard to their patch data types, patch data ghost cell widths,
    * operator stencils, etc.
    *
    * Two RefineClasses objects are consistent if they have the same number of
    * equivalence classes and each corresponding equivalence class has the same
    * characteristics as follows:
    *
    * <ul>
    *    <li> Each corresponding patch data component (d_dst, d_src, etc.)
    *         must have the same patch data type and ghost cell width.
    *    <li> d_fine_bdry_reps_var flag must have the same value.
    *    <li> The refinement operators, if any, have the same stencil width.
    *    <li> The same time interpolation operator, if any, is used.
    *    <li> The same variable fill pattern is used.
    * </ul>
    * 
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return true if test_classes is consistent with this object.
    *
    * @param[in] test_classes  RefineClasses object to check for consistency
    * @param[in] descriptor
    */
   bool
   checkConsistency(
      tbox::Pointer<RefineClasses> test_classes,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL)) const;

   /*!
    * @brief Get the size that has been allocated for the array storing refine
    * items.
    *
    * Note that this is not necessarily the same as the number of registered
    * refine items, which can be retrieved using getNumberOfRefineItems().
    * The refine item array is allocated to a default size and grown when
    * necessary or when increaseRefineItemArraySize() is called.
    */
   int
   getRefineItemArraySize() const;

   /*!
    * @brief Increase the allocated size of the array storing refine items.
    *
    * This should be used in cases where there is a large number of refine
    * items being registered with the RefineAlgorithm, to avoid frequent
    * resizing of the array.  If the size argument is less than the current
    * allocated size of the array, then the size of the array is not changed.
    *
    * @param[in] size
    */
   void
   increaseRefineItemArraySize(
      const int size);

   /*!
    * @brief Print data for all refine items to the specified output stream.
    *
    * @param[out] stream
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Print single refine item to the specified output stream.
    *
    * @param[out] stream
    * @param[in] data
    */
   void
   printRefineItem(
      std::ostream& stream,
      const RefineClasses::Data& data) const;

private:
   RefineClasses(
      const RefineClasses&);            // not implemented
   void
   operator = (
      const RefineClasses&);                     // not implemented

   /*!
    * @brief Function to compare two patch data components (with given
    * descriptor indices) for consistency.
    *
    * Two components are consistent if the are of the same patch data type and
    * have the same ghost width.
    *
    * @return true if consistent; false otherwise.
    *
    * @param[in] item_id1
    * @param[in] item_id2
    * @param[in] pd  descriptor
    */
   bool
   checkPatchDataItemConsistency(
      int item_id1,
      int item_id2,
      tbox::Pointer<hier::PatchDescriptor> pd) const;

   /*!
    * @brief Function to determine the equivalence class where a refine data
    * item belongs.
    *
    * The refine data item is compared to existing equivalence classes to
    * determine if it can be a member of any of them.
    *
    * @return If the item matches an existing equivalence class the integer
    * identifier for that equivalence class is returned.  Otherwise -1 is
    * returned.
    *
    * @param[in] data
    * @param[in] descriptor 
    */
   int
   getEquivalenceClassIndex(
      const RefineClasses::Data& data,
      tbox::Pointer<hier::PatchDescriptor> descriptor) const;

   /*!
    * The default length of the refine item array.
    */
   static int s_default_refine_item_array_size;

   /*!
    * The array of refine items.
    */
   tbox::Array<Data> d_refine_classes_data_items;

   /*!
    * The array managing equivalence classes.  Each element of the array
    * represents one equivalence class.  Each List holds integers identifying
    * which items are part of an equivalence class.  The integers index into
    * the array d_refine_classes_data_items.
    */
   tbox::Array<tbox::List<int> > d_equivalence_class_ids;

   /*!
    * The number of refine items that have been registered.
    */
   int d_num_refine_items;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/xfer/RefineClasses.I"
#endif

#endif
