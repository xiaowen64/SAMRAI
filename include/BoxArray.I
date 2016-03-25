//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BoxArray.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2142 $
// Modified:	$LastChangedDate: 2008-04-23 08:45:19 -0700 (Wed, 23 Apr 2008) $
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

#if (INCLUDE_DEPRECATED <= 2)
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
#endif

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

#if (INCLUDE_DEPRECATED <= 0)
template<int DIM> inline
Box<DIM>& BoxArray<DIM>::operator()(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_boxes.getSize()));
#endif
   return(d_boxes[i]);
}

template<int DIM> inline
const Box<DIM>& BoxArray<DIM>::operator()(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_boxes.getSize()));
#endif
   return(d_boxes[i]);
}
#endif

template <int DIM> inline
Box<DIM>& BoxArray<DIM>::operator[](const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_boxes.getSize()));
#endif
   return(d_boxes[i]);
}

template <int DIM> inline
const Box<DIM>& BoxArray<DIM>::operator[](const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_boxes.getSize()));
#endif
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

