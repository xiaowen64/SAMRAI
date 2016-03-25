//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outeredge/OuteredgeData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outeredge centered patch data type
//



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
                                         int face_normal) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   return(d_data[axis][face_normal][0].isInitialized());
}

template <int DIM, class TYPE>
inline
TYPE* OuteredgeData<DIM,TYPE>::getPointer(int axis,
                                          int face_normal,
                                          int side,
                                          int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis][face_normal][side].getPointer(depth));
}

template <int DIM, class TYPE>
inline
const TYPE* OuteredgeData<DIM,TYPE>::getPointer(int axis, 
                                                int face_normal, 
                                                int side,  
                                                int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis][face_normal][side].getPointer(depth));
}

template <int DIM, class TYPE>
inline
ArrayData<DIM,TYPE>& 
OuteredgeData<DIM,TYPE>::getArrayData(int axis,
                                      int face_normal,
                                      int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[axis][face_normal][side]);
}

template <int DIM, class TYPE>
inline
const ArrayData<DIM,TYPE>&
OuteredgeData<DIM,TYPE>::getArrayData(int axis,
                                      int face_normal,
                                      int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[axis][face_normal][side]);
}

template <int DIM, class TYPE>
inline
TYPE& OuteredgeData<DIM,TYPE>::operator()(
   const EdgeIndex<DIM>& i, 
   int depth)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif

   for (int face_normal = 0; face_normal < DIM; ++face_normal) {

      if ( face_normal != axis ) {

         for (int side = 0; side < 2; ++side) {

            if ( d_data[axis][face_normal][side].getBox().contains(i) ) {
               return ( d_data[axis][face_normal][side](i, depth));
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
              " edge index " << i << " with axis = " << axis << std::endl); 
   return ( d_data[0][0][0](i,depth) );
}

template <int DIM, class TYPE>
inline
const TYPE& OuteredgeData<DIM,TYPE>::operator()(
   const EdgeIndex<DIM>& i, 
   int depth) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif

   for (int face_normal = 0; face_normal < DIM; ++face_normal) {

      if ( face_normal != axis ) {

         for (int side = 0; side < 2; ++side) {

            if ( d_data[axis][face_normal][side].getBox().contains(i) ) {
               return ( d_data[axis][face_normal][side](i, depth));
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
              " edge index " << i << " with axis = " << axis << std::endl);
   return ( d_data[0][0][0](i,depth) );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

