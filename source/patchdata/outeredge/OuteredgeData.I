//
// File:	OuteredgeData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 692 $
// Modified:	$Date: 2005-10-28 13:37:46 -0700 (Fri, 28 Oct 2005) $
// Description:	Templated outeredge centered patch data type
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
int OuteredgeData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template <int DIM, class TYPE> inline
bool OuteredgeData<DIM,TYPE>::dataExists(int axis, 
                                         int face_nrml) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
#endif
   return(d_data[axis][face_nrml][0].isInitialized());
}

template <int DIM, class TYPE>
inline
TYPE *OuteredgeData<DIM,TYPE>::getPointer(int axis,
                                          int face_nrml,
                                          int s,
                                          int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face_nrml][s].getPointer(d));
}

template <int DIM, class TYPE>
inline
const TYPE *OuteredgeData<DIM,TYPE>::getPointer(int axis, 
                                                int face_nrml, 
                                                int s,  
                                                int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face_nrml][s].getPointer(d));
}

template <int DIM, class TYPE>
inline
ArrayData<DIM,TYPE>
&OuteredgeData<DIM,TYPE>::getArrayData(int axis,
                                       int face_nrml,
                                       int s)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
#endif
   return(d_data[axis][face_nrml][s]);
}

template <int DIM, class TYPE>
inline
const ArrayData<DIM,TYPE>
&OuteredgeData<DIM,TYPE>::getArrayData(int axis,
                                       int face_nrml,
                                       int s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
#endif
   return(d_data[axis][face_nrml][s]);
}

template <int DIM, class TYPE>
inline
TYPE& OuteredgeData<DIM,TYPE>::operator()(
   const EdgeIndex<DIM>& i, 
   int axis,
   int depth)
{
   for (int s = DIM-1; s >= 0; s--) {
      if (i[s] == d_data[axis][s][0].getBox().lower()[s] ) {
         return (d_data[axis][s][0](i,depth));
      }
      if (i[s] == d_data[axis][s][1].getBox().upper()[s] ) {
         return (d_data[axis][s][1](i,depth));
      }
   }

   /*
    * The following lines should only be executed if there's a bug
    * in the Outeredge datatype.
    */
   TBOX_ERROR("Bad index used to access outeredge data\n"
              <<"Given index is not an outeredge of this instance.\n");
   return d_data[0][0][0](i,depth);
}

template <int DIM, class TYPE>
inline
const TYPE& OuteredgeData<DIM,TYPE>::operator()(
   const EdgeIndex<DIM>& i, 
   int axis,
   int depth) const
{
   for (int s = DIM-1; s >= 0; s--) {
      if (i[s] == d_data[axis][s][0].getBox().lower()[s] ) {
         return (d_data[axis][s][0](i,depth));
      }
      if (i[s] == d_data[axis][s][1].getBox().upper()[s] ) {
         return (d_data[axis][s][1](i,depth));
      }
   }

   /*
    * The following lines should only be executed if there's a bug
    * in the Outeredge datatype.
    */
   TBOX_ERROR("Bad index used to access outeredge data\n"
              <<"Given index is not an outeredge of this instance.\n");
   return d_data[0][0][0](i,depth);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

