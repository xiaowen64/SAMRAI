//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/index/IndexData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 2251 $
// Modified:	$LastChangedDate: 2008-07-03 09:39:02 -0700 (Thu, 03 Jul 2008) $
// Description: hier::Patch data structure for irregular grid data
//

#ifndef included_pdat_IndexData
#define included_pdat_IndexData

#include "SAMRAI_config.h"

#include <vector>

#include "Box.h"
#include "PatchData.h"
#include "Index.h"
#include "IntVector.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"


namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE, class BOX_GEOMETRY> class IndexDataNode;
template<int DIM, class TYPE, class BOX_GEOMETRY> class IndexIterator;
template<int DIM, class TYPE, class BOX_GEOMETRY> class ConstIndexIterator;

/**
 * IndexData is used for storing sparse data.  The iteration over the 
 * data preserves insertion order similar to a linked list, hence data
 * may be inserted at the front (AddItem) or tail (AppendItem).
 *
 * For example, this class is used to represent embedded * boundary
 * features as a regular patch data type using the BoundaryCell class
 * as the template type.  The iteration ordering property is used to
 * visit the boundary cells in a well defined way.
 *
 * Additional insertion methods are provided for optimization in
 * specific use cases.  An extra constructor call can be avoided if
 * the "Pointer" versions of Add/AppendItem are used.   These methods
 * assume that the items are being "given" to IndexData thus IndexData
 * will delete them when needed.  The items must be dynamically 
 * allocated for this to work correctly, E.G.
 *
 * index_data.addItemPointer(someIndex, new ItemToInsert());
 *
 * The "replace" versions of the insert routines may be used when one
 * wishes to preserve the existing iteration ordering.  If an item
 * already exists the replace insert will delete the old item and put
 * the item provided in it's place (inserting it into that location in
 * the linked list).
 *
 * If one does not care about the iteration order and is frequently
 * updating the items stored, the "replace" versions are significantly
 * faster.
 *
 * The template parameter TYPE * defines the storage at each index
 * location.  IndexDaga is derived from * hier::PatchData.  
 *
 * The data type TYPE must define the following five methods which are 
 * require by this class:
 * 
 *    - \b - Default constructor (taking no arguments).  
 *    - \b - Assignment operator; i.e.,  TYPE\& operator=(const TYPE\& rhs)
 *    - \b - Copy; copySourceItem(const hier::Index<DIM>\& index,
 *                                const hier::IntVector<DIM>\& src_offset,
 *                                const TYPE\& src_item)
 *    - \b - Return size of data; size_t getDataStreamSize()
 *    - \b - Pack data into message stream; i.e., 
 *             packStream(AbstractStream\& stream,
 *    - \b - Unpack data from message stream; i.e.,
 *             unpackStream(AbstractStream\& stream,
 *             const hier::IntVector<DIM>\& offset)
 *    - \b - Write to restart;
 *             putToDatabase(tbox::Pointer<tbox::Database>\& database)
 *    - \b - Retrieve from restart;
 *             getFromDatabase(tbox::Pointer<tbox::Database>\& database)
 * 
 * The BOX_GEOMETRY template parameter defines the geometry.   BOX_GEOMETRY must
 * have a nested class name Overlap that implements he following methods:
 *
 *   - \b - getSourceOffset
 *   - \b - getDestinationBoxList
 *
 *
 * More information about the templated TYPE is provided in the IndexData
 * README file.
 * 
 * IndexData<DIM> objects are created by the IndexDataFactory<DIM>
 * factory object just as all other patch data types.
 *
 * @see pdat::IndexData
 * @see hier::PatchData
 * @see pdat::IndexDataFactory
 */

template<int DIM, class TYPE, class BOX_GEOMETRY>
class IndexData : public hier::PatchData<DIM>
{
public:
   /**
    * Define the iterator.
    */
   typedef IndexIterator<DIM,TYPE,BOX_GEOMETRY> Iterator;
   typedef ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY> ConstIterator;

   /**
    * The constructor for an IndexData object.  The box describes the interior
    * of the index space and the ghosts vector describes the ghost nodes in
    * each coordinate direction.
    */
   IndexData(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor for an IndexData object.
    */
   virtual ~IndexData<DIM,TYPE,BOX_GEOMETRY>();

   /**
    * A fast copy between the source and destination.  All data is copied
    * from the source into the destination where there is overlap in the
    * index space.
    */
   virtual void copy(const hier::PatchData<DIM>& src);
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description should have been computed
    * previously from computeIntersection().
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);
   virtual void copy2(hier::PatchData<DIM>& dst,
                      const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * Determines whether the hier::PatchData subclass can estinate the necessary
    * stream size using only index space information.
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * Pack data lying on the specified index set into the output stream.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * Unpack data from the message stream into the specified index set.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /**
    * Add a new item to the tail of the irregular index set.
    */
   void appendItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Add a pointer to a new item to the tail of the irregular index
    * set.  IndexData will delete the item when it is no longer needed
    * by IndexData.  Due to this behavior item must be dynamically
    * created (e.g. new) so that it may be deleted.
    *
    * NOTE: This is an optimization to avoid an extra constructor
    * call.  It should be used with caution, the caller MUST NOT
    * delete the referenced item.  Think of this as giving up control
    * of the item to IndexData.
    */
   void appendItemPointer(const hier::Index<DIM>& index, TYPE* item);

   /**
    * Add a new item to the head of the irregular index set
    */
   void addItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Add a pointer to a new item to the head of the irregular index
    * set.  IndexData will delete the item when it is no longer needed
    * by IndexData.  Due to this behavior item must be dynamically
    * created (e.g. new) so that it may be deleted.
    *
    * NOTE: This is an optimization to avoid an extra constructor
    * call.  It should be used with caution, the caller MUST NOT
    * delete the referenced item.  Think of this as giving up control
    * of the item to IndexData.
    */

   void addItemPointer(const hier::Index<DIM>& index, TYPE* item);


   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the orginal insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to addItem.
    */ 
   void replaceAddItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the original insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to addItemPointer.
    *
    * See addItemPointer for additional comments on pointer semantics.
    */ 
   void replaceAddItemPointer(const hier::Index<DIM>& index, TYPE* item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the orginal insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to appendItem.
    */ 
   void replaceAppendItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the original insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to appendItemPointer.
    *
    * See addItemPointer for additional comments on pointer semantics.
    */ 
   void replaceAppendItemPointer(const hier::Index<DIM>& index, TYPE* item);


   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index.
    *
    */ 
   void removeItem(const hier::Index<DIM>& index);

   /**
    * Return the number of data items (i.e. the number of indices) in   
    * the index data list.
    */ 
   int getNumberOfItems() const;

   /**
    * Remove (deallocate) any items in the irregular index set located in
    * the index space of the hier::Box.
    */
   void removeInsideBox(const hier::Box<DIM>& box);

   /**
    * Remove (deallocate) any items in the irregular index set located
    * outside of the index space of the hier::Box.
    */
   void removeOutsideBox(const hier::Box<DIM>& box);

   /**
    * Remove (deallocate) the items in the irregular index set located in 
    * the ghost region of the patch.
    */
   void removeGhostItems();

   /**
    * Remove (deallocate) all items in the irregular index set.
    */
   void removeAllItems();

   /**
    * Returns true if there is an element of the irregular index set at
    * the specified hier::Index.
    */
   bool isElement(const hier::Index<DIM>& index) const;

   /**
    * Given an index, return a pointer to the item located at that index.
    * If there is no item at the index, null is returned.
    */
   TYPE* getItem(const hier::Index<DIM>& index) const;

   /**
    * Check to make sure that the class version number is the same
    * as the restart file version number.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
      tbox::Pointer<tbox::Database> database);

   /**
    * Write out the class version number to the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase(
      tbox::Pointer<tbox::Database> database);


private:
   friend class IndexIterator<DIM,TYPE,BOX_GEOMETRY>;
   friend class ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>;

   /**
    * Returns true if element exists at offset
    */
   bool isElement(int offset) const;

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the offset.
    *
    * NOTE: This is for optimization to avoid computing 
    * offset repeatedly.
    */ 
   void removeItem(const int offset);

   /**
    * Internal routine to append item to the linked list
    * representation.
    *
    * NOTE: Offset is not strictly necessary but was include to avoid
    * computing it repeatedly.
    */
   void addItemToList(const hier::Index<DIM>& index, const int offset, TYPE& item);
   void appendItemToList(const hier::Index<DIM>& index, const int offset, TYPE& item);

   /**
    * Remove the specified node from the linked list.
    */
   void removeNodeFromList(IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *node);

   IndexData(const IndexData<DIM,TYPE,BOX_GEOMETRY>&); // not implemented
   void operator=(const IndexData<DIM,TYPE,BOX_GEOMETRY>&);	  // not implemented

   std::vector<IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *> d_data;

//   tbox::Array< IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *> d_data;

   /* 
    * Doublely linked list of nodes.
    */
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_list_head;
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_list_tail;
   int d_number_items;
};


/**
 * Class IndexDataNode holds items for the linked list.  This class should
 * be defined inside the IndexData class, but nested template classes can
 * cause some compilers to barf.  List nodes should never be seen by the
 * user of a list.
 *
 */
template<int DIM, class TYPE, class BOX_GEOMETRY>
class IndexDataNode {
public:

   friend class IndexData<DIM,TYPE,BOX_GEOMETRY>;
   friend class IndexIterator<DIM,TYPE,BOX_GEOMETRY>;
   friend class ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>;

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>();

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>(
      const hier::Index<DIM>& index,
      const int d_offset,
      TYPE& t,
      IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *n,
      IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *p);

   virtual ~IndexDataNode<DIM,TYPE,BOX_GEOMETRY>();

private:
   hier::Index<DIM> d_index;
   int   d_offset;
   TYPE* d_item;

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_next;
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_prev;
};

/**
 * Class IndexIterator is the iterator associated with the IndexData
 * This class provides methods for stepping through the
 * list that contains the irregular index set.  The user should
 * access this class through the name IndexData<DIM,TYPE>::Iterator.
 *
 * This iterator should be used as follows:
   \verbatim
   IndexData<DIM,TYPE> data;
   ...
   for (IndexData<DIM,TYPE>::Iterator iter(data); iter; iter++ {
      ... = iter();
   }
   \endverbatim
 *
 * @see tbox::List
 * @see pdat::IndexData
 * @see pdat::IndexIterator
 */

template<int DIM, class TYPE, class BOX_GEOMETRY>
class IndexIterator
{
public:
   /**
    * Default constructor for the index list iterator.  The iterator must
    * be initialized before it can be used to iterate over an IndexData object.
    */
   IndexIterator();

   /**
    * Constructor for the index list iterator.  The iterator will iterate
    * over the irregular index set of the argument data object.
    */
   IndexIterator(const IndexData<DIM,TYPE,BOX_GEOMETRY>& data);

   /**
    * Copy constructor for the index list iterator.
    */
   IndexIterator(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);

   /**
    * Assignment operator for the index list iterator.
    */
   IndexIterator<DIM,TYPE,BOX_GEOMETRY>& 
      operator=(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);

   /**
    * Destructor for the index list iterator.
    */
   ~IndexIterator<DIM,TYPE,BOX_GEOMETRY>();

   /**
    * Return the current item in the irregular index set.
    */
   TYPE& operator*();

   /**
    * Return a const reference to the current item in the irregular
    * index set.
    */
   const TYPE& operator*() const;

   /**
    * Return the current item in the irregular index set.
    */
   TYPE& operator()();

   /**
    * Return a const refeferene to the current item in the irregular
    * index set.
    */
   const TYPE& operator()() const;

   /**
    * Return the current item in the irregular index set.
    */
   TYPE& getItem();

   /**
    * Return the index of the current item in the irregular index set
    */
   const hier::Index<DIM>& getIndex() const;

   /**
    * Return true if the iterator points to a valid item in the index set.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-null if the iterator points to a valid item in the index
    * set.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid item in the index set.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next item in the index set.
    */
   void operator++(int);

   /**
    * Decrement the iterator to point to the previous item in the index set.
    */
   void operator--(int);

   /**
    * Test two iterators for equality (pointing to the same item).
    */
   bool operator==(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different items).
    */
   bool operator!=(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator) const;

   void rewindIterator();

   void fastforwardIterator();

private:
   friend class ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>;
   friend class IndexData<DIM,TYPE,BOX_GEOMETRY>;

   IndexIterator(IndexData<DIM,TYPE,BOX_GEOMETRY>* index_data,
		 IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *node);

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> &getNode();

   IndexData<DIM,TYPE,BOX_GEOMETRY> *d_index_data;

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_node;
};

#if 0

template<int DIM, class TYPE, class BOX_GEOMETRY>
class ConstIndexIterator
{
public:
   /**
    * Default constructor for the index list iterator.  The iterator must
    * be initialized before it can be used to iterate over an IndexData object.
    */
   ConstIndexIterator();

   /**
    * Constructor for the index list iterator.  The iterator will iterate
    * over the irregular index set of the argument data object.
    */
   ConstIndexIterator(const IndexData<DIM,TYPE,BOX_GEOMETRY>& data);

   /**
    * Copy constructor for the index list iterator.
    */
   ConstIndexIterator(const ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);
   ConstIndexIterator(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);

   /**
    * Assignment operator for the index list iterator.
    */
   ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& 
      operator=(const ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);
   ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& 
      operator=(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator);

   /**
    * Destructor for the index list iterator.
    */
   ~ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>();

   /**
    * Return the current item in the irregular index set.
    */
   const TYPE& operator*();

   /**
    * Return the current item in the irregular index set.
    */
   const TYPE& operator()();

   /**
    * Return the current item in the irregular index set.
    */
   const TYPE& getItem();

   /**
    * Return the index of the current item in the irregular index set
    */
   const hier::Index<DIM>& getIndex() const;

   /**
    * Return true if the iterator points to a valid item in the index set.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-null if the iterator points to a valid item in the index
    * set.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid item in the index set.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next item in the index set.
    */
   void operator++(int);

   /**
    * Test two iterators for equality (pointing to the same item).
    */
   bool operator==(const ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different items).
    */
   bool operator!=(const ConstIndexIterator<DIM,TYPE,BOX_GEOMETRY>& iterator) const;

private:
   friend class IndexIterator<DIM,TYPE,BOX_GEOMETRY>;

   const IndexData<DIM,TYPE,BOX_GEOMETRY>* d_index_data;

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *d_node;
}

#endif

    }
}

#ifndef DEBUG_NO_INLINE
#include "IndexData.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "IndexData.C"
#endif

#endif
