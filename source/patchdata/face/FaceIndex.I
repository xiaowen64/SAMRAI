//
// File:	FaceIndex.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Index for face centered patch data types
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
FaceIndex<DIM>::FaceIndex()
{
}

template<int DIM> inline
FaceIndex<DIM>::FaceIndex(const hier::Index<DIM>& rhs,
                                 const int axis,
                                 const int face)
{
   d_axis = axis;
   (*this)(0) = rhs(d_axis)+face;
   for (int i = 1; i < DIM; i++) {
      (*this)(i) = rhs((d_axis+i) % DIM);
   }
}

template<int DIM> inline
FaceIndex<DIM>::FaceIndex(const FaceIndex<DIM>& rhs)
:  hier::Index<DIM>(rhs), d_axis(rhs.d_axis)
{
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator=(const FaceIndex<DIM>& rhs)
{
   hier::Index<DIM>::operator=(rhs);
   d_axis   = rhs.d_axis;
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM>::~FaceIndex()
{
}

template<int DIM> inline
int FaceIndex<DIM>::getAxis() const
{
   return(d_axis);
}

template<int DIM> inline
void FaceIndex<DIM>::setAxis(const int axis)
{
   d_axis = axis;
}

template<int DIM> inline
hier::Index<DIM> FaceIndex<DIM>::toCell(const int face) const
{
   hier::Index<DIM> index;
   index(d_axis) = (*this)(0)+face-1;
   for (int i = 1; i < DIM; i++) {
      index((d_axis+i) % DIM) = (*this)(i);
   }
   return(index);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator+=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator+(const hier::IntVector<DIM>& rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator+=(const int rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator+(const int rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator-=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator-(const hier::IntVector<DIM>& rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator-=(const int rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator-(const int rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator*=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator*(const hier::IntVector<DIM>& rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
FaceIndex<DIM>& FaceIndex<DIM>::operator*=(const int rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
FaceIndex<DIM> FaceIndex<DIM>::operator*(const int rhs) const
{
   FaceIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
bool FaceIndex<DIM>::operator==(const FaceIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator==(rhs) && (d_axis == rhs.d_axis));
}

template<int DIM> inline
bool FaceIndex<DIM>::operator!=(const FaceIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator!=(rhs) || (d_axis != rhs.d_axis));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
