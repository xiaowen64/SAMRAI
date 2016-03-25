//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BoxArray.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	An array of boxes that complements BoxList
//

#ifndef included_hier_BoxArray_C
#define included_hier_BoxArray_C

#include "BoxArray.h"
#include "BoxList.h"
#include "Index.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "BoxArray.I"
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM>  BoxArray<DIM>::BoxArray(const int n) : d_boxes(n)
{
}

template<int DIM>  BoxArray<DIM>::BoxArray(const tbox::Array< Box<DIM> >& array)
:  d_boxes(array.getSize())
{
   const int n = array.getSize();
   for (int i = 0; i < n; i++) {
      d_boxes[i] = array[i];
   }
}

template<int DIM>  BoxArray<DIM>::BoxArray(const BoxArray<DIM>& array)
:  d_boxes(array.getNumberOfBoxes())
{
   const int n = array.getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      d_boxes[i] = array.d_boxes[i];
   }
}

template<int DIM>  BoxArray<DIM>::BoxArray(const BoxList<DIM>& list)
:  d_boxes(list.getNumberOfItems())
{
   int index = 0;
   for (typename tbox::List< Box<DIM> >::Iterator box(list); box; box++) {
      d_boxes[index++] = box();
   }
}

template<int DIM>  BoxArray<DIM>::BoxArray(const tbox::Array<tbox::DatabaseBox>& array)
{
   set_BoxArray_from_Array(array);
}

template<int DIM> BoxArray<DIM>::operator tbox::Array<tbox::DatabaseBox>() const
{
   int number_boxes = getNumberOfBoxes();
   tbox::Array<tbox::DatabaseBox> new_Array(number_boxes);

   for (int j = 0; j < number_boxes; j++) {
         new_Array[j] = (tbox::DatabaseBox) (d_boxes[j]);
   }

   return new_Array;
}

template<int DIM> void BoxArray<DIM>::set_BoxArray_from_Array(
   tbox::Array<tbox::DatabaseBox> array)
{
   const int n = array.getSize();
 
   d_boxes.resizeArray(n);
 
   for (int j = 0; j < n; j++) {
      d_boxes[j] = array[j];
   }
}

template<int DIM> BoxArray<DIM>& BoxArray<DIM>::operator=(const BoxArray<DIM>& array)
{
   if (this != &array) {
      const int n = array.getNumberOfBoxes();
      d_boxes = tbox::Array< Box<DIM> >(n);
      for (int i = 0; i < n; i++) {
         d_boxes[i] = array.d_boxes[i];
      }
   }
   return(*this);
}

template<int DIM> BoxArray<DIM>& BoxArray<DIM>::operator=(const BoxList<DIM>& list)
{
   const int n = list.getNumberOfItems();
   d_boxes = tbox::Array< Box<DIM> >(n);
   int index = 0;
   for (typename tbox::List< Box<DIM> >::Iterator box(list); box; box++) {
      d_boxes[index++] = box();
   }
   return(*this);
}

template<int DIM> bool BoxArray<DIM>::contains(const Index<DIM>& p) const
{
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      if (d_boxes[i].contains(p)) return(true);
   }
   return(false);
}

template<int DIM> void BoxArray<DIM>::grow(const IntVector<DIM>& ghosts)
{
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      d_boxes[i].grow(ghosts);
   }
}

template<int DIM> void BoxArray<DIM>::shift(const IntVector<DIM>& offset)
{
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      d_boxes[i].shift(offset);
   }
}

template<int DIM> void BoxArray<DIM>::rotate(int rotation_number)
{
   if (DIM == 2 || DIM == 3) {
      const int n = getNumberOfBoxes();
      for (int b = 0; b < n; b++) {
	 d_boxes[b].rotate(rotation_number);
      }
   } else { 
      NULL_USE(rotation_number);
      
      TBOX_ERROR("BoxArray<DIM>::rotate() error ..."
		 << "\n   Rotation only implemented for 2D and 3D " << std::endl);
   }
}

template<int DIM> void BoxArray<DIM>::refine(const IntVector<DIM>& ratio)
{
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      d_boxes[i].refine(ratio);
   }
}

template<int DIM> void BoxArray<DIM>::coarsen(const IntVector<DIM>& ratio)
{
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      d_boxes[i].coarsen(ratio);
   }
}

/*
*************************************************************************
*                                                                       *
* Return the bounding box for all boxes in the box array.               *
*                                                                       *
*************************************************************************
*/

template<int DIM> Box<DIM> BoxArray<DIM>::getBoundingBox() const
{
   Box<DIM> bbox;
   const int n = getNumberOfBoxes();
   for (int i = 0; i < n; i++) {
      bbox += d_boxes[i];
   }
   return(bbox);
}

/*
*************************************************************************
*                                                                       *
* Print boxes in array.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxArray<DIM>::print(std::ostream& os) const
{
   for (int i = 0; i < getNumberOfBoxes(); i++) {
      os << "Box # " << i << ":  " << d_boxes[i] << std::endl;
   }
}


}
}

#endif
