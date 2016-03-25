//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outernode/OuternodeData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1973 $
// Modified:	$LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
// Description:	Templated outernode centered patch data type
//


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
bool OuternodeData<DIM,TYPE>::dataExists(int face_normal) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   return(d_data[face_normal][0].isInitialized());
}

template <int DIM, class TYPE> inline
TYPE* OuternodeData<DIM,TYPE>::getPointer(int face_normal,
                                          int side,
                                          int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[face_normal][side].getPointer(depth));
}

template <int DIM, class TYPE> inline
const TYPE* OuternodeData<DIM,TYPE>::getPointer(int face_normal,
                                                int side,
                                                int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[face_normal][side].getPointer(depth));
}

template <int DIM, class TYPE> inline
ArrayData<DIM,TYPE>
&OuternodeData<DIM,TYPE>::getArrayData(int face_normal,
                                       int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[face_normal][side]);
}

template <int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>
&OuternodeData<DIM,TYPE>::getArrayData(int face_normal,
                                       int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[face_normal][side]);
}

template <int DIM, class TYPE> inline
TYPE& OuternodeData<DIM,TYPE>::operator()(
   const NodeIndex<DIM>& i, int depth)
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
   const NodeIndex<DIM>& i, int depth) const
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

