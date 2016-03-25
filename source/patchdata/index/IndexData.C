//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/index/IndexData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: hier::Patch data structure for irregular grid data
//

#ifndef included_pdat_IndexData_C
#define included_pdat_IndexData_C

#include "IndexData.h"

#include "Box.h"
#include "CellOverlap.h"
#include "tbox/Utilities.h"
#include "tbox/IOStream.h"


#define PDAT_INDEXDATA_VERSION 1

#define INDEX_NAME_BUF_SIZE 32

#ifdef DEBUG_NO_INLINE
#include "IndexData.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* The constructor for the irregular grid object simply initializes the	*
* irregular data list to be null (this is done implicitly).		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
IndexData<DIM,TYPE>::IndexData(const hier::Box<DIM>& box,
                                       const hier::IntVector<DIM>& ghosts)
:  hier::PatchData<DIM>(box, ghosts),
   d_data(hier::PatchData<DIM>::getGhostBox().size())
{

}

template<int DIM, class TYPE>
IndexData<DIM,TYPE>::~IndexData()
{

   removeAllItems();
}

/*
*************************************************************************
*									*
* The following are private and cannot be used, but they are defined	*
* here for compilers that require that every template declaration have	*
* a definition (a stupid requirement, if you ask me).			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
IndexData<DIM,TYPE>::IndexData(const IndexData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{

   // private and not used (but included for stupid compilers)
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::operator=(const IndexData<DIM,TYPE>& foo)
{
   // private and not used (but included for stupid compilers)
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Copy into dst where src overlaps on interiors.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const IndexData<DIM,TYPE> *t_src =
      dynamic_cast<const IndexData<DIM,TYPE> *>(&src);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_src != NULL);
#endif

   const hier::Box<DIM>& src_ghost_box = t_src->getGhostBox();
   removeInsideBox(src_ghost_box);
   for (typename tbox::List< hier::Index<DIM> >::Iterator
        s(t_src->d_list.listStart()); 
        s; s++) {
      appendItem(s(), *(t_src->d_data[src_ghost_box.offset(s())].d_item));
   }
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   dst.copy(*this);
}

/*
*************************************************************************
*									*
* Copy data from the source into the destination according to the	*
* overlap descriptor.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                 const hier::BoxOverlap<DIM>& overlap)
{
   const IndexData<DIM,TYPE> *t_src =
      dynamic_cast<const IndexData<DIM,TYPE> *>(&src);
   const CellOverlap<DIM> *t_overlap =
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_src != NULL);
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList();
   const hier::Box<DIM>& src_ghost_box = t_src->getGhostBox();

   for (typename hier::BoxList<DIM>::Iterator b(box_list); b; b++) {
      const hier::Box<DIM> src_box(hier::Box<DIM>::shift(b(), -src_offset));
      removeInsideBox(src_box);
      for (typename tbox::List< hier::Index<DIM> >::Iterator
           s(t_src->d_list.listStart());
           s; s++) {
         if (src_box.contains(s())) {
            TYPE new_item;
            new_item.copySourceItem(
               s(),
               src_offset,
               *(t_src->d_data[src_ghost_box.offset(s())].d_item));
            appendItem(s()+src_offset, new_item);
         }
      }
   }
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                  const hier::BoxOverlap<DIM>& overlap) const
{
   dst.copy(*this, overlap);
}

/*
*************************************************************************
*									*
* Calculate the buffer space needed to pack/unpack messages on the box	*
* region using the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool IndexData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(false);
}

template<int DIM, class TYPE>
int IndexData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const CellOverlap<DIM> *t_overlap =
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   size_t bytes = 0;
   int num_items = 0;
   const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList();
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      hier::Box<DIM> box = hier::PatchData<DIM>::getBox() *
                      hier::Box<DIM>::shift(b(), -(t_overlap->getSourceOffset()));
      for (typename hier::Box<DIM>::Iterator index(box); index; index++) {
         TYPE* item = getItem(index());
         if (item) {
            num_items++;
            bytes += item->getDataStreamSize();
         }
      }
   }
   const int index_size = DIM * tbox::AbstractStream::sizeofInt();
   bytes += (num_items * index_size + tbox::AbstractStream::sizeofInt());
   return(bytes);
}

/*
*************************************************************************
*									*
* Pack/unpack data into/out of the message streams using the index	*
* space in the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const CellOverlap<DIM>  *t_overlap =
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList();
   int num_items = 0;
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      hier::Box<DIM> box = hier::PatchData<DIM>::getBox() *
                      hier::Box<DIM>::shift(b(), -(t_overlap->getSourceOffset()));
      for (typename tbox::List< hier::Index<DIM> >::Iterator s(d_list); s; s++) {
         if (box.contains(s())) {
            num_items++;
         }
      }
   }

   stream << num_items;

   for (typename hier::BoxList<DIM>::Iterator c(boxes); c; c++) {
      hier::Box<DIM> box = hier::PatchData<DIM>::getBox() *
                      hier::Box<DIM>::shift(c(), -(t_overlap->getSourceOffset()));
      for (typename tbox::List< hier::Index<DIM> >::Iterator t(d_list); t; t++) {
         if (box.contains(t())) {
            TYPE* item = getItem(t());
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(item != NULL);
#endif
            int index_buf[DIM];
            for (int i=0; i<DIM; i++) {
               index_buf[i] = t()(i);
            }
            stream.pack(index_buf, DIM);
            item->packStream(stream);
         }
      }
   }

}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const CellOverlap<DIM> *t_overlap =
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   int num_items;
   stream >> num_items;

   const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList();
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      removeInsideBox(b());
   }

   int i;
   TYPE* items = new TYPE[num_items];
   for (i=0; i<num_items; i++) {
      int index_buf[DIM];
      stream.unpack(index_buf, DIM);
      hier::Index<DIM> index; 
      for (int j=0; j<DIM; j++) {
         index(j) = index_buf[j];
      }
      (items+i)->unpackStream(stream, t_overlap->getSourceOffset());
      addItem(index+(t_overlap->getSourceOffset()), items[i]);
   }
   delete[] items;
}

/*
*************************************************************************
*									*
* tbox::List manipulation stuff.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::appendItem(const hier::Index<DIM>& index,
                                       const TYPE& item)
{
   if (hier::PatchData<DIM>::getGhostBox().contains(index)) {
      if (isElement(index)) {
         removeItem(index);
      }
      TYPE* new_item = new TYPE();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(new_item != NULL);
#endif
      *new_item = item; 
      d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item =  new_item;
      d_list.appendItem(index);
   }
}
 

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::addItem(const hier::Index<DIM>& index, const TYPE& item)
{
   if (hier::PatchData<DIM>::getGhostBox().contains(index)) {
      if (isElement(index)) {
         removeItem(index);
      }
      TYPE* new_item = new TYPE();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(new_item != NULL);
#endif
      *new_item = item;
      d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item =  new_item;
      d_list.addItem(index);
   }
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::removeItem(const hier::Index<DIM>& index)
{
   bool found = false;
   for (typename tbox::List< hier::Index<DIM> >::Iterator l(d_list.listStart());
        (!found) && l; l++) {
      if (l() == index) {
         d_list.removeItem(l);
         found = true;
      }
   }
   if (found) {
      TYPE* tmp = d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item;
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(tmp != NULL);
#endif
      d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item = NULL;
      delete tmp;
   } 
}

template <int DIM, class TYPE>
int IndexData<DIM,TYPE>::getNumberItems() const
{
   return(d_list.getNumberItems());
}


template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::removeInsideBox(const hier::Box<DIM>& box)
{
   typename tbox::List< hier::Index<DIM> >::Iterator l(d_list);
   while (l) {
      if (box.contains(l())) {
         hier::Index<DIM> index(l());
         l++;
         removeItem(index);
      } else {
         l++;
      }
   }
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::removeOutsideBox(const hier::Box<DIM>& box)
{
   typename tbox::List< hier::Index<DIM> >::Iterator l(d_list);
   while (l) {
      if (!box.contains(l())) {
         hier::Index<DIM> index(l());
         l++;
         removeItem(index);;
      } else {
         l++;
      }
   }
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::removeGhostItems()
{
   removeOutsideBox(hier::PatchData<DIM>::getBox());
}

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::removeAllItems()
{
   removeInsideBox(hier::PatchData<DIM>::getGhostBox());
}

template<int DIM, class TYPE>
bool IndexData<DIM,TYPE>::isElement(const hier::Index<DIM>& index) const
{
   bool is_element = false;
   if ((hier::PatchData<DIM>::getGhostBox().contains(index)) && 
      (d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item != NULL)) {
      is_element = true;
   }
   return(is_element);
}

/*
*************************************************************************
*									*
* Just checks to make sure that the class version is the same		*
* as the restart file version number.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_INDEXDATA_VERSION");
   if (ver != PDAT_INDEXDATA_VERSION){
      TBOX_ERROR("IndexData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl); 
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( INDEX_NAME_BUF_SIZE > (5 + 1 + 6 + 1) );
#endif

   char index_keyword[INDEX_NAME_BUF_SIZE];

   int item_count = 0;
   bool item_found = true;

   do {
      sprintf(index_keyword,"index_data_%06d",item_count);

      if (database->isDatabase(index_keyword)) {

         tbox::Pointer<tbox::Database> item_db =
            database->getDatabase(index_keyword);

         tbox::Array<int> index_array = item_db->getIntegerArray(index_keyword);
         hier::Index<DIM> index;
         for (int j = 0; j < DIM; j++) {
            index(j) = index_array[j];
         }

         TYPE item;
         item.getFromDatabase(item_db);

         appendItem(index, item);

      } else {
         item_found = false;
      }

      item_count++;

   } while (item_found);

}

/*
*************************************************************************
*									*
* Just writes out the class version number to the database.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void IndexData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_INDEXDATA_VERSION",PDAT_INDEXDATA_VERSION);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( INDEX_NAME_BUF_SIZE > (5 + 1 + 6 + 1) );
#endif

   char index_keyword[INDEX_NAME_BUF_SIZE]; 
   int item_count = 0;
   for (typename tbox::List< hier::Index<DIM> >::Iterator s(d_list.listStart()); s; s++) {
      sprintf(index_keyword,"index_data_%06d",item_count);
      hier::Index<DIM> index = s();
      tbox::Array<int> index_array(DIM);
      for (int i = 0; i < DIM; i++) {
         index_array[i] = index(i);
      }

      tbox::Pointer<tbox::Database> item_db =
         database->putDatabase(index_keyword);

      item_db->putIntegerArray(index_keyword, index_array);

      TYPE* item = getItem(index);

      item->putToDatabase(item_db);

      item_count++;
   }
}

template<int DIM, class TYPE>
TYPE* IndexData<DIM,TYPE>::getItem(const hier::Index<DIM>& index) const
{
   TYPE* item;
   if (!isElement(index)) {
      item = NULL;
   } else {
      item = d_data[hier::PatchData<DIM>::getGhostBox().offset(index)].d_item;
   }

   return item;
}

template<int DIM, class TYPE>
IndexIterator<DIM,TYPE>::IndexIterator()
{
}

template<int DIM, class TYPE>
IndexIterator<DIM,TYPE>::IndexIterator(
   IndexData<DIM,TYPE>& data)
{
   d_index_data = &data;
   d_iterator = data.d_list;
}

template<int DIM, class TYPE>
IndexIterator<DIM,TYPE>::IndexIterator(
   const IndexIterator<DIM,TYPE>& iter)
:  d_iterator(iter.d_iterator),
   d_index_data(iter.d_index_data)
{
}

#ifdef ENABLE_CONST_ITERATOR
// ConstIndexIterator
template<int DIM, class TYPE>
ConstIndexIterator<DIM,TYPE>::ConstIndexIterator()
{
}

template<int DIM, class TYPE>
ConstIndexIterator<DIM,TYPE>::ConstIndexIterator(
   const IndexData<DIM,TYPE>& data)
{
   d_index_data = &data;
   d_iterator = data.d_list;
}

template<int DIM, class TYPE>
ConstIndexIterator<DIM,TYPE>::ConstIndexIterator(
   const ConstIndexIterator<DIM,TYPE>& iter)
:  d_iterator(iter.d_iterator),
   d_index_data(iter.d_index_data)
{
}

template<int DIM, class TYPE>
ConstIndexIterator<DIM,TYPE>::ConstIndexIterator(
   const IndexIterator<DIM,TYPE>& iter)
:  d_iterator(iter.d_iterator),
   d_index_data(iter.d_index_data)
{
}
#endif


template<int DIM, class TYPE>
IndexDataNode<DIM,TYPE>::IndexDataNode()
{
   d_item = NULL;
}

template<int DIM, class TYPE>
IndexDataNode<DIM,TYPE>::~IndexDataNode()
{
}


}
}

#endif
