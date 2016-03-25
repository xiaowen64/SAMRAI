//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/cell/CellData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated cell centered patch data type
//


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
TYPE* CellData<DIM,TYPE>::getPointer(int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data.getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* CellData<DIM,TYPE>::getPointer(int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data.getPointer(depth));
}

template<int DIM, class TYPE> inline
TYPE& CellData<DIM,TYPE>::operator()(const CellIndex<DIM>& i, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data(i,depth));
}

template<int DIM, class TYPE> inline
const TYPE& CellData<DIM,TYPE>::operator()(const CellIndex<DIM>& i, int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data(i,depth));
}

template<int DIM, class TYPE> inline
void CellData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, d);
}

template<int DIM, class TYPE> inline
void CellData<DIM,TYPE>::fill(const TYPE& t,
                              const hier::Box<DIM>& box,
                              int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, box, d);
}

template<int DIM, class TYPE> inline 
void CellData<DIM,TYPE>::fillAll(const TYPE& t)
{
   d_data.fillAll(t);
}

template<int DIM, class TYPE>
inline void CellData<DIM,TYPE>::fillAll(const TYPE& t, 
                                        const hier::Box<DIM>& box)
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
