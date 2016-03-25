//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/index/IndexData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: hier::Patch data structure for irregular grid data
//

#ifndef included_pdat_IndexData
#define included_pdat_IndexData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_hier_Index
#include "Index.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#define ENABLE_CONST_ITERATOR 1

namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> class IndexDataNode;
template<int DIM, class TYPE> class IndexIterator;
#ifdef ENABLE_CONST_ITERATOR
template<int DIM, class TYPE> class ConstIndexIterator;
#endif

/**
 * Class IndexData<DIM> is a templated patch data type for manipulating 
 * data over irregular cell-centered index sets.  The template parameter TYPE 
 * defines the storage at each index location.  It is derived from 
 * hier::PatchData<DIM>.   For example, this class is used to represent embedded 
 * boundary features as a regular patch data type using the BoundaryCell class 
 * as the template type.
 *
 * The data type TYPE must define the following five methods which are 
 * require by this class:
 * 


 *    - \b - Default constructor (taking no arguments).  
 *    - \b - Assignment operator; i.e.,  TYPE\& operator=(const TYPE\& rhs)
 *    - \b - Copy; copySourceItem(const hier::Index<DIM>\& index,
                                    const hier::IntVector<DIM>\& src_offset,
                                    const TYPE\& src_item)
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

template<int DIM, class TYPE>
class IndexData : public hier::PatchData<DIM>
{
public:
   /**
    * Define the iterator.
    */
   typedef IndexIterator<DIM,TYPE> Iterator;
#ifdef ENABLE_CONST_ITERATOR
   typedef ConstIndexIterator<DIM,TYPE> ConstIterator;
#endif

   /**
    * The constructor for an IndexData object.  The box describes the interior
    * of the index space and the ghosts vector describes the ghost nodes in
    * each coordinate direction.
    */
   IndexData(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor for an IndexData object.
    */
   virtual ~IndexData<DIM,TYPE>();

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
    * Add a new item to the tail of the irregular index set
    */
   void appendItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Add a new item to the head of the irregular index set
    */
   void addItem(const hier::Index<DIM>& index, const TYPE& item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index.
    */ 
   void removeItem(const hier::Index<DIM>& index);

   /**
    * Return the number of data items (i.e. the number of indices) in   
    * the index data list.
    */ 
   int getNumberItems() const;

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

   friend class IndexIterator<DIM,TYPE>;
#ifdef ENABLE_CONST_ITERATOR
   friend class ConstIndexIterator<DIM,TYPE>;
#endif

private:
   IndexData(const IndexData<DIM,TYPE>&); // not implemented
   void operator=(const IndexData<DIM,TYPE>&);	  // not implemented

   tbox::List< hier::Index<DIM> > d_list;
   tbox::Array< IndexDataNode<DIM,TYPE> > d_data;
};

template<int DIM, class TYPE>
class IndexDataNode {
public:

   friend class IndexData<DIM,TYPE>;
   friend class IndexIterator<DIM,TYPE>;
#ifdef ENABLE_CONST_ITERATOR
   friend class ConstIndexIterator<DIM,TYPE>;
#endif

   IndexDataNode();

   virtual ~IndexDataNode<DIM,TYPE>();

private:
   TYPE* d_item;
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

template<int DIM, class TYPE>
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
   IndexIterator(IndexData<DIM,TYPE>& data);

   /**
    * Copy constructor for the index list iterator.
    */
   IndexIterator(const IndexIterator<DIM,TYPE>& iterator);

   /**
    * Assignment operator for the index list iterator.
    */
   IndexIterator<DIM,TYPE>& operator=(const IndexIterator<DIM,TYPE>& iterator);

   /**
    * Destructor for the index list iterator.
    */
   ~IndexIterator<DIM,TYPE>();

   /**
    * Return the current item in the irregular index set.
    */
   TYPE& operator*();

   /**
    * Return the current item in the irregular index set.
    */
   TYPE& operator()();

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
    * Test two iterators for equality (pointing to the same item).
    */
   bool operator==(const IndexIterator<DIM,TYPE>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different items).
    */
   bool operator!=(const IndexIterator<DIM,TYPE>& iterator) const;

#ifdef ENABLE_CONST_ITERATOR
   friend class ConstIndexIterator<DIM,TYPE>;
#endif
private:
   typename tbox::List< hier::Index<DIM> >::Iterator d_iterator;
   IndexData<DIM,TYPE>* d_index_data;
};

#ifdef ENABLE_CONST_ITERATOR
template<int DIM, class TYPE>
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
   ConstIndexIterator(const IndexData<DIM,TYPE>& data);

   /**
    * Copy constructor for the index list iterator.
    */
   ConstIndexIterator(const ConstIndexIterator<DIM,TYPE>& iterator);
   ConstIndexIterator(const IndexIterator<DIM,TYPE>& iterator);

   /**
    * Assignment operator for the index list iterator.
    */
   ConstIndexIterator<DIM,TYPE>& operator=(const ConstIndexIterator<DIM,TYPE>& iterator);
   ConstIndexIterator<DIM,TYPE>& operator=(const IndexIterator<DIM,TYPE>& iterator);

   /**
    * Destructor for the index list iterator.
    */
   ~ConstIndexIterator<DIM,TYPE>();

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
   bool operator==(const ConstIndexIterator<DIM,TYPE>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different items).
    */
   bool operator!=(const ConstIndexIterator<DIM,TYPE>& iterator) const;

   friend class IndexIterator<DIM,TYPE>;
private:
   typename tbox::List< hier::Index<DIM> >::Iterator d_iterator;
   const IndexData<DIM,TYPE>* d_index_data;
};
#endif

}
}

#ifndef DEBUG_NO_INLINE
#include "IndexData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "IndexData.C"
#endif
