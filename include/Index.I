//
// File:	Index.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Interface for the AMR Index object
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> inline
Index<DIM>::Index() : IntVector<DIM>()
{
}

template<int DIM> inline
Index<DIM>::Index(const int i) : IntVector<DIM>(i)
{
}

#if INCLUDE_DEPRECATED < 2
template<int DIM> inline
Index<DIM>::Index(const int i, const int j) : IntVector<DIM>(i,j)
{
}

template<int DIM> inline
Index<DIM>::Index(const int i, const int j, const int k)
:  IntVector<DIM>(i,j,k)
{
}
#endif

template<int DIM> inline
Index<DIM>::Index(const tbox::Array<int> i)
:  IntVector<DIM>(i)
{
}

template<int DIM> inline
Index<DIM>::Index(const Index<DIM>& rhs) : IntVector<DIM>(rhs)
{
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator=(const Index<DIM>& rhs)
{
   IntVector<DIM>::operator=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM>::~Index()
{
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator+=(const IntVector<DIM>& rhs)
{
   IntVector<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator+(const IntVector<DIM>& rhs) const
{
   Index<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator+=(const int rhs)
{
   IntVector<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator+(const int rhs) const
{
   Index<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator-=(const IntVector<DIM>& rhs)
{
   IntVector<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator-(const IntVector<DIM>& rhs) const
{
   Index<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator-=(const int rhs)
{
   IntVector<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator-(const int rhs) const
{
   Index<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator*=(const IntVector<DIM>& rhs)
{
   IntVector<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator*(const IntVector<DIM>& rhs) const
{
   Index<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator*=(const int rhs)
{
   IntVector<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator*(const int rhs) const
{
   Index<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator/=(const IntVector<DIM>& rhs)
{
   IntVector<DIM>::operator/=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator/(const IntVector<DIM>& rhs) const
{
   Index<DIM> tmp = *this;
   tmp /= rhs;
   return(tmp);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator/=(const int rhs)
{
   IntVector<DIM>::operator/=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM> Index<DIM>::operator/(const int rhs) const
{
   Index<DIM> tmp = *this;
   tmp /= rhs;
   return(tmp);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

