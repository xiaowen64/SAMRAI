//
// File:	OuternodeData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated outernode centered patch data type
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

template <int DIM, class TYPE> inline
int OuternodeData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template <int DIM, class TYPE> inline
bool OuternodeData<DIM,TYPE>::dataExists(const int axis) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
#endif
   return(d_data[axis][0].isInitialized());
}

template <int DIM, class TYPE> inline
TYPE *OuternodeData<DIM,TYPE>::getPointer(const int axis,
                                          const int side,
                                          const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side].getPointer(d));
}

template <int DIM, class TYPE> inline
const TYPE *OuternodeData<DIM,TYPE>::getPointer(const int axis,
                                                  const int side,
                                                  const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side].getPointer(d));
}

template <int DIM, class TYPE> inline
ArrayData<DIM,TYPE>
&OuternodeData<DIM,TYPE>::getArrayData(const int axis,
                                         const int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
#endif
   return(d_data[axis][side]);
}

template <int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>
&OuternodeData<DIM,TYPE>::getArrayData(const int axis,
                                         const int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
#endif
   return(d_data[axis][side]);
}

template <int DIM, class TYPE> inline
TYPE& OuternodeData<DIM,TYPE>::operator()(
   const NodeIndex<DIM>& i, const int depth)
{
   for (int d = DIM-1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d] ) {
         return (d_data[d][0](i,depth));
      }
      if (i[d] == d_data[d][1].getBox().upper()[d] ) {
         return (d_data[d][1](i,depth));
      }
   }

   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
              <<"Given index is not an outernode of this instance.\n");
   return d_data[0][0](i,depth);
}

template <int DIM, class TYPE> inline
const TYPE& OuternodeData<DIM,TYPE>::operator()(
   const NodeIndex<DIM>& i, const int depth) const
{
   for (int d = DIM-1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d] ) {
         return (d_data[d][0](i,depth));
      }
      if (i[d] == d_data[d][1].getBox().upper()[d] ) {
         return (d_data[d][1](i,depth));
      }
   }
   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
              <<"Given index is not an outernode of this instance.\n");
   return d_data[0][0](i,depth);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

