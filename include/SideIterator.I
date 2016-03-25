//
// File:	SideIterator.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Iterator for side centered patch data types
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
SideIterator<DIM>::SideIterator()
{
}

template<int DIM> inline
SideIterator<DIM>::SideIterator(const hier::Box<DIM>& box, const int axis)
:  d_index(box.lower(), axis, SideIndex<DIM>::Lower),
   d_box(SideGeometry<DIM>::toSideBox(box, axis))
{
}

template<int DIM> inline
SideIterator<DIM>::SideIterator(const SideIterator<DIM>& iter)
:  d_index(iter.d_index), d_box(iter.d_box)
{
}

template<int DIM> inline
SideIterator<DIM>&
SideIterator<DIM>::operator=(const SideIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
SideIterator<DIM>::~SideIterator()
{
}

template<int DIM> inline
const SideIndex<DIM>& SideIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const SideIndex<DIM>& SideIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
SideIterator<DIM>::operator bool() const
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
SideIterator<DIM>::operator const void*() const
{
   return(SideIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool SideIterator<DIM>::operator!() const
{
   return(!SideIterator<DIM>::operator bool());
}
   
template<int DIM> inline
void SideIterator<DIM>::operator++(int)
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
bool SideIterator<DIM>::operator==(const SideIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool SideIterator<DIM>::operator!=(const SideIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
