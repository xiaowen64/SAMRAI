//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/side/SideIndex.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Index for side centered patch data types
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
SideIndex<DIM>::SideIndex()
{
}

template<int DIM> inline
SideIndex<DIM>::SideIndex(const hier::Index<DIM>& rhs,
                                 const int axis,
                                 const int side)
{
   d_axis = axis;

   for (int i = 0; i < DIM; i++) {
      (*this)(i) = rhs(i);
   }

   (*this)(d_axis) += side;
}

template<int DIM> inline
SideIndex<DIM>::SideIndex(const SideIndex<DIM>& rhs)
:  hier::Index<DIM>(rhs), d_axis(rhs.d_axis)
{
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator=(const SideIndex<DIM>& rhs)
{
   hier::Index<DIM>::operator=(rhs);
   d_axis   = rhs.d_axis;
   return(*this);
}

template<int DIM> inline
SideIndex<DIM>::~SideIndex()
{
}

template<int DIM> inline
int SideIndex<DIM>::getAxis() const
{
   return(d_axis);
}

template<int DIM> inline
void SideIndex<DIM>::setAxis(const int axis)
{
   d_axis = axis;
}

template<int DIM> inline
hier::Index<DIM> SideIndex<DIM>::toCell(const int side) const
{
   hier::Index<DIM> index;

   for (int i = 0; i < DIM; i++) {
      index(i) = (*this)(i);
   }

   index(d_axis) += (side - 1);

   return(index);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator+=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator+(const hier::IntVector<DIM>& rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator+=(const int rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator+(const int rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator-=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator-(const hier::IntVector<DIM>& rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator-=(const int rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator-(const int rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator*=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator*(const hier::IntVector<DIM>& rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
SideIndex<DIM>& SideIndex<DIM>::operator*=(const int rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
SideIndex<DIM> SideIndex<DIM>::operator*(const int rhs) const
{
   SideIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
bool SideIndex<DIM>::operator==(const SideIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator==(rhs) && (d_axis == rhs.d_axis));
}

template<int DIM> inline
bool SideIndex<DIM>::operator!=(const SideIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator!=(rhs) || (d_axis != rhs.d_axis));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
