//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/cell/CellIndex.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Index for face centered patch data types
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
CellIndex<DIM>::CellIndex()
{
}

template<int DIM> inline
CellIndex<DIM>::CellIndex(const hier::Index<DIM>& rhs) : hier::Index<DIM>(rhs)
{
}

template<int DIM> inline
CellIndex<DIM>::CellIndex(const CellIndex<DIM>& rhs) : hier::Index<DIM>(rhs)
{
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator=(const CellIndex<DIM>& rhs)
{
   hier::Index<DIM>::operator=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM>::~CellIndex()
{
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator+=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator+(const hier::IntVector<DIM>& rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator+=(const int rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator+(const int rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator-=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator-(const hier::IntVector<DIM>& rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator-=(const int rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator-(const int rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator*=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator*(const hier::IntVector<DIM>& rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
CellIndex<DIM>& CellIndex<DIM>::operator*=(const int rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
CellIndex<DIM> CellIndex<DIM>::operator*(const int rhs) const
{
   CellIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
bool CellIndex<DIM>::operator==(const CellIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator==(rhs));
}

template<int DIM> inline
bool CellIndex<DIM>::operator!=(const CellIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator!=(rhs));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
