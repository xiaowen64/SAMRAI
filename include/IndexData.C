//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/index/IndexData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description: hier::Patch data structure for irregular grid data
//

#ifndef included_pdat_IndexData_C
#define included_pdat_IndexData_C

#include "IndexData.h"

#include "Box.h"
#include "BoxOverlap.h"
#include "tbox/Utilities.h"
#include "tbox/IOStream.h"

#define PDAT_INDEXDATA_VERSION 1

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

template<int DIM, class TYPE, class BOX_GEOMETRY>
IndexData<DIM,TYPE,BOX_GEOMETRY>::IndexData(const hier::Box<DIM>& box,
                                       const hier::IntVector<DIM>& ghosts)
:  hier::PatchData<DIM>(box, ghosts),
   d_data(hier::PatchData<DIM>::getGhostBox().size()),
   d_list_head(NULL),
   d_list_tail(NULL),
   d_number_items(0)
{

}

template<int DIM, class TYPE, class BOX_GEOMETRY>
IndexData<DIM,TYPE,BOX_GEOMETRY>::~IndexData()
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

template<int DIM, class TYPE, class BOX_GEOMETRY>
IndexData<DIM,TYPE,BOX_GEOMETRY>::IndexData(const IndexData<DIM,TYPE,BOX_GEOMETRY>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{

   // private and not used (but included for some compilers)
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::operator=(const IndexData<DIM,TYPE,BOX_GEOMETRY>& foo)
{
   // private and not used (but included for some compilers)
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Copy into dst where src overlaps on interiors.			*
*									*
*************************************************************************
*/
template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::copy(const hier::PatchData<DIM>& src)
{
   const IndexData<DIM,TYPE,BOX_GEOMETRY> *t_src =
      dynamic_cast<const IndexData<DIM,TYPE,BOX_GEOMETRY> *>(&src);

   TBOX_CHECK_ASSERT(t_src != NULL);

   const hier::Box<DIM>& src_ghost_box = t_src->getGhostBox();
   removeInsideBox(src_ghost_box);
   for (typename IndexData<DIM, TYPE, BOX_GEOMETRY>::Iterator 
	   s(*t_src); 
        s; 
	s++) {
      if (this -> getGhostBox().contains(s.getNode().d_index)) {
	 appendItem(s.getNode().d_index, *s.getNode().d_item);
      }
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::copy2(hier::PatchData<DIM>& dst) const
{
   dst.copy(*this);
}

/*
*************************************************************************
*									*
* Copy data from the source into the destination according to the	*
* overlap descriptor.							*
*                                                                       *
*************************************************************************
*/
 	
template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::copy(const hier::PatchData<DIM>& src,
                               const hier::BoxOverlap<DIM>& overlap)
{
   const IndexData<DIM,TYPE,BOX_GEOMETRY> *t_src =
      dynamic_cast<const IndexData<DIM,TYPE,BOX_GEOMETRY> *>(&src);
   const typename BOX_GEOMETRY::Overlap *t_overlap = 
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);

   TBOX_CHECK_ASSERT(t_src != NULL);
   TBOX_CHECK_ASSERT(t_overlap != NULL);

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList();
   const hier::Box<DIM>& src_ghost_box = t_src->getGhostBox();

   for (typename hier::BoxList<DIM>::Iterator b(box_list); b; b++) {
      const hier::Box<DIM>& dst_box = b();
      const hier::Box<DIM> src_box(hier::Box<DIM>::shift(b(), -src_offset));
      removeInsideBox(dst_box);
      for (typename IndexData<DIM, TYPE, BOX_GEOMETRY>::Iterator s(*t_src);
           s; 
	   s++) {
         if (src_box.contains(s.getNode().d_index)) {
            TYPE new_item;
            new_item.copySourceItem(
               s.getNode().d_index,
               src_offset,
               *(t_src->d_data[src_ghost_box.offset(s.getNode().d_index)] -> d_item));
            appendItem(s.getNode().d_index+src_offset, new_item);
         }
      }
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::copy2(hier::PatchData<DIM>& dst,
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

template<int DIM, class TYPE, class BOX_GEOMETRY>
bool IndexData<DIM,TYPE,BOX_GEOMETRY>::canEstimateStreamSizeFromBox() const
{
   return(false);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
int IndexData<DIM,TYPE,BOX_GEOMETRY>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const typename BOX_GEOMETRY::Overlap *t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_CHECK_ASSERT(t_overlap != NULL);

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

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const typename BOX_GEOMETRY::Overlap *t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_CHECK_ASSERT(t_overlap != NULL);

   const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList();
   int num_items = 0;
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      hier::Box<DIM> box = hier::PatchData<DIM>::getBox() *
                      hier::Box<DIM>::shift(b(), -(t_overlap->getSourceOffset()));
      for (typename IndexData<DIM,TYPE,BOX_GEOMETRY>::Iterator s(*this); s; s++) {
         if (box.contains(s.getNode().d_index)) {
            num_items++;
         }
      }
   }

   stream << num_items;

   for (typename hier::BoxList<DIM>::Iterator c(boxes); c; c++) {
      hier::Box<DIM> box = hier::PatchData<DIM>::getBox() *
                      hier::Box<DIM>::shift(c(), -(t_overlap->getSourceOffset()));
      for (typename IndexData<DIM,TYPE,BOX_GEOMETRY>::Iterator t(*this); t; t++) {
         if (box.contains(t.getNode().d_index)) {
            TYPE* item = &t();
            TBOX_CHECK_ASSERT(item != NULL);

            int index_buf[DIM];
            for (int i=0; i<DIM; i++) {
               index_buf[i] = t.getNode().d_index(i);
            }
            stream.pack(index_buf, DIM);
            item->packStream(stream);
         }
      }
   }

}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const typename BOX_GEOMETRY::Overlap *t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_CHECK_ASSERT(t_overlap != NULL);

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

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::appendItem(const hier::Index<DIM>& index,
                                       const TYPE& item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }

   TYPE* new_item = new TYPE();
   TBOX_CHECK_ASSERT(new_item != NULL);

   *new_item = item; 
   addItemToList(index, offset, *new_item);
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::appendItemPointer(const hier::Index<DIM>& index,
                                       TYPE* item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   appendItemToList(index, offset, *item);
}
 

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::addItem(const hier::Index<DIM>& index, const TYPE& item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   TYPE* new_item = new TYPE();
   TBOX_CHECK_ASSERT(new_item != NULL);

   *new_item = item;
   addItemToList(index, offset, *new_item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::addItemPointer(const hier::Index<DIM>& index, TYPE* item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   addItemToList(index, offset, *item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::replaceAddItem(const hier::Index<DIM>& index, const  TYPE& item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_CHECK_ASSERT(new_item != NULL);

   *new_item = item;

   if(node == NULL) {

      addItemToList(index, offset, *new_item);

   } else {
      delete node -> d_item;

      node -> d_item = new_item;
   }
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::replaceAddItemPointer(const hier::Index<DIM>& index, TYPE* item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>* node = d_data[offset];

   if(node == NULL) {

      addItemToList(index, offset, *item);

   } else {

      delete node -> d_item;

      node -> d_item = item;
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::replaceAppendItem(const hier::Index<DIM>& index, const  TYPE& item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_CHECK_ASSERT(new_item != NULL);

   *new_item = item;

   if(node == NULL) {

      appendItemToList(index, offset, *new_item);

   } else {
      delete node -> d_item;

      node -> d_item = new_item;
   }
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::replaceAppendItemPointer(const hier::Index<DIM>& index, TYPE* item)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>* node = d_data[offset];

   if(node == NULL) {

      appendItemToList(index, offset, *item);

   } else {

      delete node -> d_item;

      node -> d_item = item;
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeItem(const hier::Index<DIM>& index)
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));
   
   int offset = hier::PatchData<DIM>::getGhostBox().offset(index);
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   removeItem(offset);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeItem(const int offset)
{
   TBOX_CHECK_ASSERT(offset >= 0 && offset <= hier::PatchData<DIM>::getGhostBox().size());

   IndexDataNode<DIM,TYPE,BOX_GEOMETRY>* node = d_data[offset];

   TBOX_CHECK_ASSERT(node);

   removeNodeFromList(node);

   delete node -> d_item;
   delete node;

   d_data[offset] = NULL;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::addItemToList(const hier::Index<DIM>& index, const int offset, TYPE& item)
{
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *new_node = 
      new IndexDataNode<DIM,TYPE,BOX_GEOMETRY>(index, offset, item, d_list_head, NULL);

   if (d_list_head) {
      d_list_head->d_prev = new_node;
   }

   d_list_head = new_node;

   if (!d_list_tail) {
      d_list_tail = new_node;
   }

   d_data[offset] = new_node;

   d_number_items++;
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::appendItemToList(const hier::Index<DIM>& index, const int offset, TYPE& item)
{
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *new_node = 
      new IndexDataNode<DIM,TYPE,BOX_GEOMETRY>(index, offset, item, NULL, d_list_tail);

   if (d_list_tail) {
      d_list_tail->d_next = new_node;
   }

   d_list_tail = new_node;

   if (!d_list_head) {
      d_list_head = new_node;
   }

   d_data[offset] = new_node;

   d_number_items++;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeNodeFromList(IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *node)
{
   if ((d_list_head == node) && (d_list_tail == node)) {
      d_list_head = d_list_tail = NULL;
      
   } else if (d_list_head == node) {
      d_list_head = node->d_next;
      node->d_next->d_prev = NULL;

   } else if (d_list_tail == node) {
      d_list_tail = node->d_prev;
      node->d_prev->d_next = NULL;

   } else {
      node->d_next->d_prev = node->d_prev;
      node->d_prev->d_next = node->d_next;
   }

   d_data[node -> d_offset] = NULL;

   d_number_items--;
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
int IndexData<DIM,TYPE,BOX_GEOMETRY>::getNumberOfItems() const
{
   return(d_number_items);
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeInsideBox(const hier::Box<DIM>& box)
{
   typename IndexData<DIM, TYPE, BOX_GEOMETRY>::Iterator l(*this);
   
   while (l) {
      if (box.contains(l.getNode().d_index)) {
         hier::Index<DIM> index(l.getNode().d_index);
         l++;
         removeItem(index);
      } else {
         l++;
      }
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeOutsideBox(const hier::Box<DIM>& box)
{
   typename IndexData<DIM, TYPE, BOX_GEOMETRY>::Iterator l(*this);

   while (l) {
      if (!box.contains(l.getNode().d_index)) {
         hier::Index<DIM> index(l.getNode().d_index);
         l++;
         removeItem(index);;
      } else {
         l++;
      }
   }
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeGhostItems()
{
   removeOutsideBox(hier::PatchData<DIM>::getBox());
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::removeAllItems()
{
   removeInsideBox(hier::PatchData<DIM>::getGhostBox());
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
bool IndexData<DIM,TYPE,BOX_GEOMETRY>::isElement(const hier::Index<DIM>& index) const
{
   TBOX_CHECK_ASSERT(hier::PatchData<DIM>::getGhostBox().contains(index));

   return d_data[hier::PatchData<DIM>::getGhostBox().offset(index)] != NULL;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
bool IndexData<DIM,TYPE,BOX_GEOMETRY>::isElement(int offset) const
{
   return d_data[offset] != NULL;
}


/*
*************************************************************************
*									*
* Just checks to make sure that the class version is the same		*
* as the restart file version number.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
   TBOX_CHECK_ASSERT(!database.isNull());

   int ver = database->getInteger("PDAT_INDEXDATA_VERSION");
   if (ver != PDAT_INDEXDATA_VERSION){
      TBOX_ERROR("IndexData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl); 
   }


   int item_count = 0;
   bool item_found = true;

   do {
      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(item_count, 6);

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

template<int DIM, class TYPE, class BOX_GEOMETRY>
void IndexData<DIM,TYPE,BOX_GEOMETRY>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
   TBOX_CHECK_ASSERT(!database.isNull());

   database->putInteger("PDAT_INDEXDATA_VERSION",PDAT_INDEXDATA_VERSION);

   int item_count = 0;
   for (typename IndexData<DIM, TYPE, BOX_GEOMETRY>::Iterator s(*this); s; s++) {
      
      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(item_count, 6);
      hier::Index<DIM> index = s.getNode().d_index;
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

template<int DIM, class TYPE, class BOX_GEOMETRY>
TYPE* IndexData<DIM,TYPE,BOX_GEOMETRY>::getItem(const hier::Index<DIM>& index) const
{
   TYPE* item;
   if (!isElement(index)) {
      item = NULL;
   } else {
      item = d_data[hier::PatchData<DIM>::getGhostBox().offset(index)] -> d_item;
   }

   return item;
}

}
}

#endif
