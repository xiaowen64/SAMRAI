/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes. 
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenClasses
#define included_xfer_CoarsenClasses

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <iostream>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Maintain a collection of coarsen items and organize them
 * into equivalence classes.
 * 
 * CoarsenClasses is used by the CoarsenSchedule and CoarsenAlgorithm
 * lasses to manage a set of coarsening data items that describe
 * coarsening of patch data between two levels in an AMR hierarchy.
 * Specifically, this class organizes these items into equivalence clases, so
 * that items are grouped together if they are considered equivalent.
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
 *   <li> They have the same VariableFillPattern
 * </ul>
 */

class CoarsenClasses:public tbox::DescribedClass
{
public:
   /*!
    * @brief Nested class used to describe a coarsening operation
    * between patch data components on an AMR hierarchy.
    */
   class Data
   {
public:
      /*!
       * @brief Destination patch data component
       */
      int d_dst;

      /*!
       * @brief Source patch data component
       */
      int d_src;

      /*!
       * @brief Boolean flag that is set to true when it is desired that fine
       * data values have priority over coarse data values when data exists
       * at the same location in the mesh on levels with different resolutions.
       */
      bool d_fine_bdry_reps_var;

      /*!
       * If the coarsen operation requires data to be coarsened from the fine
       * level's ghost regions onto the coarse data representing the same
       * mesh space, then this IntVector tells how wide the ghost region to
       * coarsen will be.  It is represented in terms of the coarse level's
       * index space.
       */
      hier::IntVector d_gcw_to_coarsen;

      /*!
       * @brief Coarsening operator.
       */
      tbox::Pointer<hier::CoarsenOperator> d_opcoarsen;

      /*!
       * @brief Identifier of equivalence class where this item belongs.  All
       * items of the same equivalence class will have the same value.
       */
      int d_class_id;

      /*!
       * @brief An array index telling where this item sits in an array of
       * coarsen items.
       */
      int d_tag;

      /*!
       * @brief VariableFillPattern that can restrict the stencil of the data
       * coarsened by the CoarsenSchedule.
       */
      tbox::Pointer<VariableFillPattern> d_var_fill_pattern;

      /*!
       * @brief Constructor.
       *
       * @param[in] dim Dimension.
       */
      explicit Data(
         tbox::Dimension dim):
         d_gcw_to_coarsen(dim) {
      }

private:
      Data();  //not implemented 
   };

   /*!
    * @brief The constructor creates an empty array of coarsen classes.
    *
    * @deprecated fill_coarse_data is no longer used.
    */
   explicit CoarsenClasses(
      bool fill_coarse_data);

   /*!
    * @brief The virtual destructor destroys the coarsen data items owned
    * by this object.
    */
   virtual ~CoarsenClasses();

   /*!
    * Return number of equivalence classes maintained by this object.
    */
   int
   getNumberOfEquivalenceClasses() const;

   /*!
    * @brief Return total number of coarsen items that have been registered
    * and stored in the CoarsenClasses object
    */
   int
   getNumberOfCoarsenItems() const;

   /*!
    * @brief Return number of coarsen data items in the equivalence class
    * represented by the given integer identifier.
    *
    * @param[in] equiv_class_id
    */
   int
   getNumberOfItemsInEquivalenceClass(
      int equiv_class_id) const;

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
   const CoarsenClasses::Data&
   getClassRepresentative(
      int equiv_class_id) const;

   /*!
    * @brief Get a coarsen item from the array of all coarsen items held by
    * this object.
    *
    * The internal storage of the coarsen items held by this class is not
    * controlled by the user, so this method is intended for use when looping
    * over all of the items, from 0 to getNumberOfCoarsenItems()-1, or when
    * looping over the integers in the List obtained from getIterator().
    *
    * @return A coarsen item identified by an integer id.
    *
    * @param[in] coarsen_item_array_id
    */
   CoarsenClasses::Data&
   getCoarsenItem(
      const int coarsen_item_array_id);

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
    * ids that can be passed into getCoarsenItem().  The array ids in a
    * single list all correspond to coarsen items in a single equivalence
    * class.
    *
    * @param[in] equiv_class_id
    */
   tbox::List<int>::Iterator
   getIterator(
      int equiv_class_id);

   /*!
    * @brief Give a data item to the CoarsenClasses object, which will store
    * it with the proper equivalence class.
    *
    * If the item belongs in an existing equivalence class, it will be added
    * there, otherwise a new equivalence class will be created for this item.
    * The internal data of the item will be changed so that it stores an
    * integer identifier of its equivalence class.
    *
    * An error will occur with a descriptive message if the data item is
    * not valid.  See checkCoarsenItem() for explanation of validity.
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
      CoarsenClasses::Data& data,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL));

   /*!
    * Check coarsen data item for validity.
    *
    * A coarsen data item is invalid if any of its patch data components are
    * negative, or if its source data does not have sufficient ghost width
    * for the stencil of the coarsen operator. or if the data types of
    * the source and destination data are not compatible to be copied
    * from one to another.
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
   checkCoarsenItem(
      const CoarsenClasses::Data& data_item,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL)) const;

   /*!
    * @brief Compare CoarsenClasses object with another CoarsenClasses object.
    *
    * This method checks if the equivalence classes held by the two objects
    * match with regard to their patch data types, patch data ghost cell widths,    * operator stencils, etc.
    *
    * Two CoarsenClasses objects are consistent if they have the same number of
    * equivalence classes and each corresponding equivalence class has the same
    * characteristics as follows:
    *
    * <ul>
    *    <li> Each corresponding patch data component (d_dst and d_src)
    *         must have the same patch data type and ghost cell width.
    *    <li> d_fine_bdry_reps_var flag must have the same value.
    *    <li> The coarsening operators, if any, have the same stencil width.
    *    <li> The same variable fill pattern is used.
    * </ul>
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return true if test_classes is consistent with this object.
    *
    * @param[in] test_classes  CoarsenClasses object to check for consistency
    * @param[in] descriptor
    */
   bool
   checkConsistency(
      tbox::Pointer<CoarsenClasses> test_classes,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL)) const;

   /*!
    * @brief Get the size that has been allocated for the array storing coarsen
    * items.
    *
    * Note that this is not necessarily the same as the number of registered
    * coarsen items, which can be retrieved using getNumberOfCoarsenItems().
    * The coarsen item array is allocated to a default size and grown when
    * necessary or when increaseCoarsenItemArraySize() is called.
    */
   int
   getCoarsenItemArraySize() const;

   /*!
    * @brief Increase the allocated size of the array storing coarsen items.
    *
    * This should be used in cases where there is a large number of coarsen
    * items being registered with the CoarsenAlgorithm, to avoid frequent
    * resizing of the array.  If the size argument is less than the current
    * allocated size of the array, then the size of the array is not changed.
    *
    * @param[in] size
    */
   void
   increaseCoarsenItemArraySize(
      const int size,
      const tbox::Dimension& dim);

   /*!
    * @brief Print data for all coarsen items to the specified output stream.
    *
    * @param[out] stream
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Print single coarsen item to the specified output stream.
    *
    * @param[out] stream
    * @param[in] data
    */
   void
   printCoarsenItem(
      std::ostream& stream,
      const CoarsenClasses::Data& data) const;

private:
   CoarsenClasses(
      const CoarsenClasses&);                   // not implemented
   void
   operator = (
      const CoarsenClasses&);                     // not implemented

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
    * @brief Function to determine the equivalence class where a coarsen data
    * item belongs.
    *
    * The coarsen data item is compared to existing equivalence classes to
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
      const CoarsenClasses::Data& data,
      tbox::Pointer<hier::PatchDescriptor> descriptor =
         tbox::Pointer<hier::PatchDescriptor>(NULL)) const;

   /*!
    * The default length of the coarsen item array.
    */
   static int s_default_coarsen_item_array_size;

   /*!
    * @deprecated  No longer used
    */
   bool d_fill_coarse_data;

   /*!
    * The array of coarsen items.
    */
   tbox::Array<CoarsenClasses::Data> d_coarsen_classes_data_items;

   /*!
    * The array managing equivalence classes.  Each element of the array
    * represents one equivalence class.  Each List holds integers identifying
    * which items are part of an equivalence class.  The integers index into
    * the array d_coarsen_classes_data_items.
    */
   tbox::Array<tbox::List<int> > d_equivalence_class_ids;

   /*!
    * The number of coarsen items that have been registered.
    */
   int d_num_coarsen_items;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/xfer/CoarsenClasses.I"
#endif

#endif
