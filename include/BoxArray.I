//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxArray.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	An array of boxes that complements BoxList
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> inline
BoxArray<DIM>::~BoxArray()
{
}

template<int DIM> inline
int BoxArray<DIM>::getNumberOfBoxes() const
{
   return(d_boxes.getSize());
}

template<int DIM> inline
int BoxArray<DIM>::size() const
{
   return(d_boxes.getSize());
}

template<int DIM> inline
Box<DIM>& BoxArray<DIM>::getBox(const int i)
{
   return(d_boxes[i]);
}

template<int DIM> inline
const Box<DIM>& BoxArray<DIM>::getBox(const int i) const
{
   return(d_boxes[i]);
}

template<int DIM> inline
BoxArray<DIM>& BoxArray<DIM>::operator=(const tbox::Array<tbox::DatabaseBox>& array)
{
   set_BoxArray_from_Array(array);
   return (*this);
}

template<int DIM> inline
BoxArray<DIM>& BoxArray<DIM>::BoxArray_from_Array(
   tbox::Array<tbox::DatabaseBox> array)
{
   set_BoxArray_from_Array(array);
   return (*this);
}

template<int DIM> inline
Box<DIM>& BoxArray<DIM>::operator()(const int i)
{
   return(d_boxes[i]);
}

template<int DIM> inline
const Box<DIM>& BoxArray<DIM>::operator()(const int i) const
{
   return(d_boxes[i]);
}

template<int DIM> inline
void BoxArray<DIM>::resizeBoxArray(const int n)
{
   d_boxes.resizeArray(n);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

