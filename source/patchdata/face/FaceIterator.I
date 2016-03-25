//
// File:	FaceIterator.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Iterator for face centered patch data types
//

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
FaceIterator<DIM>::FaceIterator()
{
}

template<int DIM> inline
FaceIterator<DIM>::FaceIterator(const hier::Box<DIM>& box, const int axis)
:  d_index(box.lower(), axis, FaceIndex<DIM>::Lower),
   d_box(FaceGeometry<DIM>::toFaceBox(box, axis))
{
}

template<int DIM> inline
FaceIterator<DIM>::FaceIterator(const FaceIterator<DIM>& iter)
:  d_index(iter.d_index), d_box(iter.d_box)
{
}

template<int DIM> inline
FaceIterator<DIM>&
FaceIterator<DIM>::operator=(const FaceIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
FaceIterator<DIM>::~FaceIterator()
{
}

template<int DIM> inline
const FaceIndex<DIM>& FaceIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const FaceIndex<DIM>& FaceIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
FaceIterator<DIM>::operator bool() const
{
   bool retval = true;
   for (int i = 0; i < DIM; i ++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }

   return(retval);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM> inline
FaceIterator<DIM>::operator const void*() const
{
   return(FaceIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool FaceIterator<DIM>::operator!() const
{
   return(!FaceIterator<DIM>::operator bool());
}
   
template<int DIM> inline
void FaceIterator<DIM>::operator++(int)
{
   d_index(0)++;
   for (int i = 0; i < DIM-1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i+1)++;
      } else {
         break;
      }
   }
}

template<int DIM> inline
bool FaceIterator<DIM>::operator==(const FaceIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool FaceIterator<DIM>::operator!=(const FaceIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
