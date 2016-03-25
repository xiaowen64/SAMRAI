//
// File:	CellData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated cell centered patch data type
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

template<int DIM, class TYPE> inline
int CellData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>& CellData<DIM,TYPE>::getArrayData()
{
   return(d_data);
}

template<int DIM, class TYPE>inline
const ArrayData<DIM,TYPE>& CellData<DIM, TYPE>::getArrayData() const
{
   return(d_data);
}

template<int DIM, class TYPE> inline
TYPE *CellData<DIM,TYPE>::getPointer(const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data.getPointer(d));
}

template<int DIM, class TYPE> inline
const TYPE *CellData<DIM,TYPE>::getPointer(const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data.getPointer(d));
}

template<int DIM, class TYPE> inline
TYPE& CellData<DIM,TYPE>::operator()(const CellIndex<DIM>& i, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data(i,d));
}

template<int DIM, class TYPE> inline
const TYPE& CellData<DIM,TYPE>::operator()(const CellIndex<DIM>& i, const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data(i,d));
}

template<int DIM, class TYPE> inline
void CellData<DIM,TYPE>::fill(const TYPE& t, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, d);
}

template<int DIM, class TYPE> inline
void CellData<DIM,TYPE>::fill(const TYPE& t,
                                const hier::Box<DIM>& box,
                                const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, box, d);
}

template<int DIM, class TYPE> inline 
void CellData<DIM,TYPE>::fillAll(const TYPE& t)
{
   d_data.fillAll(t);
}

template<int DIM, class TYPE>
inline void CellData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   d_data.fillAll(t, box);
}

template<int DIM, class TYPE>
inline void CellData<DIM,TYPE>::copyOnBox(
   const CellData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{
   const hier::Box<DIM> cell_box = CellGeometry<DIM>::toCellBox(box);
   d_data.copy(src.getArrayData(), cell_box);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif  
}
}
