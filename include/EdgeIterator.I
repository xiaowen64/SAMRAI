//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/edge/EdgeIterator.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Iterator for edge centered patch data types
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
EdgeIterator<DIM>::EdgeIterator()
{
}

template<int DIM> inline
EdgeIterator<DIM>::EdgeIterator(const hier::Box<DIM>& box, const int axis)
: d_box(EdgeGeometry<DIM>::toEdgeBox(box,axis)),
  d_index(box.lower(), axis, 0)
{
}

template<int DIM> inline
EdgeIterator<DIM>::EdgeIterator(const EdgeIterator<DIM>& iter)
:  d_box(iter.d_box), d_index(iter.d_index)
{
}

template<int DIM> inline
EdgeIterator<DIM>&
EdgeIterator<DIM>::operator=(const EdgeIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
EdgeIterator<DIM>::~EdgeIterator()
{
}

template<int DIM> inline
const EdgeIndex<DIM>& EdgeIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const EdgeIndex<DIM>& EdgeIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
EdgeIterator<DIM>::operator bool() const
{
   bool retval = true;
   for (int i = 0; i < DIM; i++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }
   return(retval);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM> inline
EdgeIterator<DIM>::operator const void*() const
{
   return(EdgeIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool EdgeIterator<DIM>::operator!() const
{
   return(!EdgeIterator<DIM>::operator bool());
}
   
template<int DIM> inline
void EdgeIterator<DIM>::operator++(int)
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
bool EdgeIterator<DIM>::operator==(const EdgeIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool EdgeIterator<DIM>::operator!=(const EdgeIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
