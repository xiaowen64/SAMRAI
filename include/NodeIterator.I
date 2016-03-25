//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeIterator.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Iterator for node centered patch data types
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
NodeIterator<DIM>::NodeIterator()
{
}

template<int DIM> inline
NodeIterator<DIM>::NodeIterator(const hier::Box<DIM>& box)
:  d_index(box.lower(), hier::IntVector<DIM>(0)),
   d_box(NodeGeometry<DIM>::toNodeBox(box))
{
}

template<int DIM> inline
NodeIterator<DIM>::NodeIterator(const NodeIterator<DIM>& iter)
:  d_index(iter.d_index), d_box(iter.d_box)
{
}

template<int DIM> inline
NodeIterator<DIM>&
NodeIterator<DIM>::operator=(const NodeIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
NodeIterator<DIM>::~NodeIterator()
{
}

template<int DIM> inline
const NodeIndex<DIM>& NodeIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const NodeIndex<DIM>& NodeIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
NodeIterator<DIM>::operator bool() const
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
NodeIterator<DIM>::operator const void*() const
{
   return(NodeIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool NodeIterator<DIM>::operator!() const
{
   return(!NodeIterator<DIM>::operator bool());
}
   
template<int DIM> inline
void NodeIterator<DIM>::operator++(int)
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
bool NodeIterator<DIM>::operator==(const NodeIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool NodeIterator<DIM>::operator!=(const NodeIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
