//
// File:	Box.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Box representing a portion of the AMR index space
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif
#ifndef included_tbox_IEEE
#include "tbox/IEEE.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
Box<DIM>::Box() : d_lo(tbox::IEEE::getINT_MAX()), d_hi(tbox::IEEE::getINT_MIN())
{
}

template<int DIM> inline
Box<DIM>::Box(const Index<DIM>& lower, const Index<DIM>& upper)
:  d_lo(lower), d_hi(upper)
{
}

template<int DIM> inline
Box<DIM>::Box(const Box<DIM>& box) : d_lo(box.d_lo), d_hi(box.d_hi)
{
}

template<int DIM> inline
Box<DIM>::Box(const tbox::DatabaseBox& box)
{
   set_Box_from_DatabaseBox(box);
}

template<int DIM> inline
Box<DIM>& Box<DIM>::Box_from_DatabaseBox(const tbox::DatabaseBox& box)
{
   set_Box_from_DatabaseBox(box);
   return (*this);
}

template<int DIM> inline
Box<DIM>& Box<DIM>::operator=(const Box<DIM>& box)
{
   d_lo = box.d_lo;
   d_hi = box.d_hi;
   return(*this);
}

template<int DIM> inline
Box<DIM>& Box<DIM>::operator=(const tbox::DatabaseBox& box)
{
   return Box_from_DatabaseBox(box);
}

template<int DIM> inline
Box<DIM>::~Box()
{
}

template<int DIM> inline
Index<DIM>& Box<DIM>::lower()
{
   return(d_lo);
}

template<int DIM> inline
Index<DIM>& Box<DIM>::upper()
{
   return(d_hi);
}

template<int DIM> inline
const Index<DIM>& Box<DIM>::lower() const
{
   return(d_lo);
}

template<int DIM> inline
const Index<DIM>& Box<DIM>::upper() const
{
   return(d_hi);
}

template<int DIM> inline
int& Box<DIM>::lower(const int i)
{
   return(d_lo(i));
}

template<int DIM> inline
int& Box<DIM>::upper(const int i)
{
   return(d_hi(i));
}

template<int DIM> inline
const int& Box<DIM>::lower(const int i) const
{
   return(d_lo(i));
}

template<int DIM> inline
const int& Box<DIM>::upper(const int i) const
{
   return(d_hi(i));
}

template<int DIM> inline
bool Box<DIM>::empty() const
{
   for (int i = 0; i < DIM; i++) {
      if (d_hi(i) < d_lo(i)) return(true);
   }

   return(false);
}

template<int DIM> inline
void Box<DIM>::setEmpty()
{
   d_lo = Index<DIM>(tbox::IEEE::getINT_MAX());
   d_hi = Index<DIM>(tbox::IEEE::getINT_MIN());
}

template<int DIM> inline 
IntVector<DIM> Box<DIM>::numberCells() const
{
   if (empty()) {
      return(IntVector<DIM>(0));
   } else {
      return(d_hi-d_lo+1);
   }
}
 
template<int DIM> inline 
int Box<DIM>::numberCells(const int i) const
{
   if (empty()) {
      return(0);
   } else {
      return(d_hi(i)-d_lo(i)+1);
   }
}

template<int DIM> inline
int Box<DIM>::size() const
{
   int mysize = 0;
   if (!empty()) {
      mysize = 1;
      for (int i = 0; i < DIM; i++) 
	 mysize *= (d_hi(i)-d_lo(i)+1);
   }

   return mysize;
}


template<int DIM> inline
int Box<DIM>::offset(const Index<DIM>& p) const
{
   int myoffset = 0;

   for (int i = DIM-1; i > 0; i--) {
      myoffset = (d_hi(i-1)-d_lo(i-1)+1)*(p(i)-d_lo(i) + myoffset);
   }
   myoffset += p(0)-d_lo(0);

   return(myoffset);
}


template<int DIM> inline
bool Box<DIM>::contains(const Index<DIM>& p) const
{
   for (int i = 0; i < DIM; i++) {
      if ( (p(i) < d_lo(i)) || (p(i) > d_hi(i)) ) return(false);
   }
   return(true);
}

template<int DIM> inline
Box<DIM>::operator tbox::DatabaseBox()
{
   return DatabaseBox_from_Box();
}

template<int DIM> inline
Box<DIM>::operator tbox::DatabaseBox() const
{
   return DatabaseBox_from_Box();
}

template<int DIM> inline
int Box<DIM>::operator==(const Box<DIM>& box) const
{
   return(((d_lo==box.d_lo) && (d_hi==box.d_hi)) || (empty() && box.empty()));
}

template<int DIM> inline
int Box<DIM>::operator!=(const Box<DIM>& box) const
{
   return(!(*this == box));
}

template<int DIM> inline
Box<DIM> Box<DIM>::operator*(const Box<DIM>& box) const
{
   Box<DIM> both(*this);
   both.d_lo.max(box.d_lo);
   both.d_hi.min(box.d_hi);
   return(both);
}

template<int DIM> inline
bool Box<DIM>::intersects(const Box<DIM>& box) const
{
   for (int i = 0; i < DIM; i++) {
      if (  tbox::Utilities::imax(d_lo(i), box.d_lo(i))
          > tbox::Utilities::imin(d_hi(i), box.d_hi(i)) ) return(false);
   }
   return(true);
}

template<int DIM> inline
Box<DIM> Box<DIM>::operator+(const Box<DIM>& box) const
{
   Box<DIM> bbox(*this);
   bbox += box;
   return(bbox);
}

template<int DIM> inline
void Box<DIM>::grow(const IntVector<DIM>& ghosts)
{
   if (!empty()) {
      d_lo -= ghosts;
      d_hi += ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::grow(const int direction, const int ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (direction >= 0) && (direction < DIM) );
#endif
   if (!empty()) {
      d_lo(direction) -= ghosts;
      d_hi(direction) += ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::growLower(const IntVector<DIM>& ghosts)
{
   if (!empty()) {
      d_lo -= ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::growLower(const int direction, const int ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (direction >= 0) && (direction < DIM) );
#endif
   if (!empty()) {
      d_lo(direction) -= ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::growUpper(const IntVector<DIM>& ghosts)
{
   if (!empty()) {
      d_hi += ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::growUpper(const int direction, const int ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (direction >= 0) && (direction < DIM) );
#endif
   if (!empty()) {
      d_hi(direction) += ghosts;
   }
}

template<int DIM> inline
void Box<DIM>::shift(const IntVector<DIM>& offset)
{
   d_lo += offset;
   d_hi += offset;
}

template<int DIM> inline
void Box<DIM>::shift(const int direction, const int offset)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (direction >= 0) && (direction < DIM) );
#endif
   d_lo(direction) += offset;
   d_hi(direction) += offset;
}

template<int DIM> inline
void Box<DIM>::refine(const IntVector<DIM>& ratio)
{
   bool negative_ratio = false;
   for (int d = 0; d < DIM; d++) {
      if (ratio(d) < 0) {
         negative_ratio = true;
         break;
      }
   } 

   if (!negative_ratio) {
      d_lo *= ratio;
      d_hi = d_hi*ratio + (ratio-1);
   } else {
      for (int i = 0; i < DIM; i++) {
         if (ratio(i) > 0) {
            d_lo(i) *= ratio(i);
            d_hi(i) = d_hi(i)*ratio(i) + (ratio(i)-1);
         } else {
            d_lo(i) = coarsen(d_lo(i), -ratio(i));
            d_hi(i) = coarsen(d_hi(i), -ratio(i));
         }
      }
   }
}

template<int DIM> inline
int Box<DIM>::coarsen(const int index, const int ratio)
{
   return(index < 0 ? (index+1)/ratio-1 : index/ratio);
}

template<int DIM> inline
void Box<DIM>::coarsen(const IntVector<DIM>& ratio)
{
   for (int i = 0; i < DIM; i++) {
      d_lo(i) = coarsen(d_lo(i), ratio(i));
      d_hi(i) = coarsen(d_hi(i), ratio(i));
   }
}

template<int DIM> inline
Box<DIM> Box<DIM>::grow(const Box<DIM>& box, const IntVector<DIM>& ghosts)
{
   Box<DIM> tmp = box;
   tmp.grow(ghosts);
   return(tmp);
}

template<int DIM> inline
Box<DIM> Box<DIM>::shift(const Box<DIM>& box, const IntVector<DIM>& ghosts)
{
   return(Box<DIM>(box.lower()+ghosts, box.upper()+ghosts));
}

template<int DIM> inline
Box<DIM> Box<DIM>::refine(const Box<DIM>& box, const IntVector<DIM>& ratio)
{
   Box<DIM> tmp = box;
   tmp.refine(ratio);
   return(tmp);
}

template<int DIM> inline
Box<DIM> Box<DIM>::coarsen(const Box<DIM>& box, const IntVector<DIM>& ratio)
{
   Box<DIM> tmp = box;
   tmp.coarsen(ratio);
   return(tmp);
}

template<int DIM> inline
BoxIterator<DIM>::BoxIterator()
{
}

template<int DIM> inline
BoxIterator<DIM>::BoxIterator(const Box<DIM>& box)
:  d_index(box.lower()), d_box(box)
{
}

template<int DIM> inline
BoxIterator<DIM>::BoxIterator(const BoxIterator<DIM>& iter)
:  d_index(iter.d_index), d_box(iter.d_box)
{
}

template<int DIM> inline
void BoxIterator<DIM>::initialize(const Box<DIM>& box)
{
   d_index = box.lower();
   d_box = box;
   return;
}

template<int DIM> inline
BoxIterator<DIM>&
BoxIterator<DIM>::operator=(const BoxIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
BoxIterator<DIM>::~BoxIterator()
{
}

template<int DIM> inline
const Index<DIM>& BoxIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const Index<DIM>& BoxIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
BoxIterator<DIM>::operator bool() const
{
   for (int i = 0; i < DIM; i++) {
      if (d_index(i) > d_box.upper(i)) return(false);
   }
   return(true);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION 
template<int DIM> inline
BoxIterator<DIM>::operator const void*() const
{
   return(BoxIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool BoxIterator<DIM>::operator!() const
{
   return(!BoxIterator<DIM>::operator bool());
}

template<int DIM> inline
void BoxIterator<DIM>::operator++(int)
{
   d_index(0)++;
   for (int i = 0; i < (DIM-1); i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i+1)++;
      } else 
         break;
   }
}

template<int DIM> inline
bool BoxIterator<DIM>::operator==(const BoxIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool BoxIterator<DIM>::operator!=(const BoxIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

