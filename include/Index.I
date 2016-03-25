//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/Index.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2147 $
// Modified:	$LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
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

template<int DIM> inline
Index<DIM>::Index(const int i, const int j) : IntVector<DIM>(i,j)
{
}

template<int DIM> inline
Index<DIM>::Index(const int i, const int j, const int k)
:  IntVector<DIM>(i,j,k)
{
}

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
Index<DIM>::Index(const IntVector<DIM>& rhs) : IntVector<DIM>(rhs)
{
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator=(const Index<DIM>& rhs)
{
   IntVector<DIM>::operator=(rhs);
   return(*this);
}

template<int DIM> inline
Index<DIM>& Index<DIM>::operator=(const IntVector<DIM>& rhs)
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

