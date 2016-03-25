//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/edge/EdgeIndex.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Index for edge centered patch data types
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
EdgeIndex<DIM>::EdgeIndex()
{
}

template<int DIM> inline
EdgeIndex<DIM>::EdgeIndex(const hier::Index<DIM>& rhs,
                                 const int axis,
                                 const int edge)
{
   d_axis = axis;
   for (int i = 0; i < DIM; i++) {
      (*this)(i) = rhs(i);
   }

   if (DIM > 1) {
      (*this)((d_axis+1) % DIM) += edge%2;
   }
   for (int j = 2; j < DIM; j++) {
      (*this)((d_axis+j) % DIM) += (edge/(1 << (j-1)))%2;
   }
}

template<int DIM> inline
EdgeIndex<DIM>::EdgeIndex(const EdgeIndex<DIM>& rhs)
:  hier::Index<DIM>(rhs), d_axis(rhs.d_axis)
{
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator=(const EdgeIndex<DIM>& rhs)
{
   hier::Index<DIM>::operator=(rhs);
   d_axis   = rhs.d_axis;
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM>::~EdgeIndex()
{
}

template<int DIM> inline
int EdgeIndex<DIM>::getAxis() const
{
   return(d_axis);
}

template<int DIM> inline
void EdgeIndex<DIM>::setAxis(const int axis)
{
   d_axis = axis;
}

template<int DIM> inline
hier::Index<DIM> EdgeIndex<DIM>::toCell(const int edge) const
{
   hier::Index<DIM> index;

   for (int i = 0; i < DIM; i++) {
      index(i) = (*this)(i);
   }

   if (DIM > 1) {
      index((d_axis+1) % DIM) += ((edge%2) - 1);
   }
   for (int j = 2; j < DIM; j++) {
      index((d_axis+j) % DIM) += (((edge/(2 << (j-1)))%2) - 1);
   }
   return(index);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator+=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator+(const hier::IntVector<DIM>& rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator+=(const int rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator+(const int rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator-=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator-(const hier::IntVector<DIM>& rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator-=(const int rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator-(const int rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator*=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator*(const hier::IntVector<DIM>& rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
EdgeIndex<DIM>& EdgeIndex<DIM>::operator*=(const int rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
EdgeIndex<DIM> EdgeIndex<DIM>::operator*(const int rhs) const
{
   EdgeIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
bool EdgeIndex<DIM>::operator==(const EdgeIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator==(rhs) && (d_axis == rhs.d_axis));
}

template<int DIM> inline
bool EdgeIndex<DIM>::operator!=(const EdgeIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator!=(rhs) || (d_axis != rhs.d_axis));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
