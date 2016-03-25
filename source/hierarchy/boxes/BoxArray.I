//
// File:	BoxArray.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

