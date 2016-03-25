//
// File:	IntVector.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	A n-dimensional integer vector
//

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include <iostream>
using namespace std;

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> inline 
IntVector<DIM>::IntVector()
{
}

template<int DIM> inline 
IntVector<DIM>::IntVector(const int value)
{
   for(int i=0; i < DIM; i++)
      d_vector[i] = value;
}


#if INCLUDE_DEPRECATED < 2
// SGS These are for backward compatiblity
// Is there a way to signal this at compile time...seem to recall
// seeing that done somewhere in STL code.
template<int DIM> inline 
IntVector<DIM>::IntVector(const int i, const int j)
{
   if (DIM != 2 ) {
      TBOX_ERROR("Error in IntVector::IntVector...\n" 
 		 << "Attempting to use 2D constructor for " 
 		 << DIM << "d  IntVector" 
 		 << endl);
   }
   else {
      // The "funny" indexing prevents a warning when compiling for 
      // DIM < 2.  This code is only reached if DIM = 2 when 
      // executing.
      d_vector[0]               = i;
      d_vector[DIM > 1 ? 1 : 0] = j;
   }
}

template<int DIM> inline
IntVector<DIM>::IntVector(const int i, const int j, const int k)
{
   if (DIM != 3 ) {
      TBOX_ERROR("Error in IntVector::IntVector...\n" 
 		 << "Attempting to use 3D constructor for " 
 		 << DIM << "d  IntVector" 
 		 << endl);
   } else {
      // The "funny" indexing prevents a warning when compiling for 
      // DIM < 2.  This code is only reached if DIM = 3 when 
      // executing.
      d_vector[0]               = i;
      d_vector[DIM > 1 ? 1 : 0] = j;
      d_vector[DIM > 2 ? 2 : 0] = k;
   }
}
#endif

template<int DIM> inline
IntVector<DIM>::IntVector(const tbox::Array<int> a)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(a.getSize() == DIM);
#endif

   for (int i = 0; i < DIM; i++) 
      d_vector[i] = a[i];
}

template<int DIM> inline
IntVector<DIM>::IntVector(const IntVector& rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] = rhs.d_vector[i];
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator=(const IntVector<DIM>& rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] = rhs.d_vector[i];

   return(*this);
}

template<int DIM> inline
IntVector<DIM>::~IntVector()
{
}

template<int DIM> inline
int& IntVector<DIM>::operator()(const int i)
{
   return(d_vector[i]);
}

template<int DIM> inline
const int& IntVector<DIM>::operator()(const int i) const
{
   return(d_vector[i]);
}

template<int DIM> inline
IntVector<DIM>::operator int*()
{
   return(d_vector);
}

template<int DIM> inline
IntVector<DIM>::operator const int*() const
{
   return(d_vector);
}


template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator+=(const IntVector<DIM>& rhs)
{
   for(int i = 0; i < DIM; i++)
	   d_vector[i] += rhs.d_vector[i];

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator+(const IntVector<DIM>& rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator+=(const int rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] += rhs;

   return(*this);
}


template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator+(const int rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator-=(const IntVector<DIM>& rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] -= rhs.d_vector[i];

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator-(const IntVector<DIM>& rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator-=(const int rhs)
{
   for(int i=0; i < DIM; i++)
      d_vector[i] -= rhs;

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator-(const int rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator*=(const IntVector<DIM>& rhs)
{
   for(int i=0; i < DIM; i++)
      d_vector[i] *= rhs.d_vector[i];

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator*(const IntVector<DIM>& rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator*=(const int rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] *= rhs;

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator*(const int rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator/=(const IntVector<DIM>& rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] /= rhs.d_vector[i];

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator/(const IntVector<DIM>& rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp /= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM>& IntVector<DIM>::operator/=(const int rhs)
{
   for(int i = 0; i < DIM; i++)
      d_vector[i] /= rhs;

   return(*this);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator/(const int rhs) const
{
   IntVector<DIM> tmp = *this;
   tmp /= rhs;
   return(tmp);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::operator-() const
{
   IntVector<DIM> tmp;
   for(int i = 0; i < DIM; i++)
      tmp.d_vector[i] = -d_vector[i];

   return(tmp);
}

template<int DIM> inline
bool IntVector<DIM>::operator==(const IntVector<DIM>& rhs) const
{
   bool result = true;

   for(int i = 0; result && (i < DIM); i++)
      result = result && (d_vector[i] == rhs.d_vector[i]);

   return result;
}

template<int DIM> inline
bool IntVector<DIM>::operator!=(const IntVector<DIM>& rhs) const
{
   return(!(*this == rhs));
}

template<int DIM> inline
void IntVector<DIM>::min(IntVector<DIM> const & rhs)
{

   for(int i = 0; i < DIM; i++)
      if (rhs.d_vector[i] < d_vector[i]) 
	 d_vector[i] = rhs.d_vector[i];

}

template<int DIM> inline
bool IntVector<DIM>::operator<(const IntVector<DIM>& rhs) const
{

   bool result = true;

   for(int i = 0; result && (i < DIM); i++)
      result = result && (d_vector[i] < rhs.d_vector[i]);

   return result;
}

template<int DIM> inline
bool IntVector<DIM>::operator<=(const IntVector<DIM>& rhs) const
{
   bool result = true;

   for(int i = 0; result && (i < DIM); i++)
      result = result && (d_vector[i] <= rhs.d_vector[i]);

   return result;
	 
}

template<int DIM> inline
bool IntVector<DIM>::operator>(const IntVector<DIM>& rhs) const
{
   bool result = true;

   for(int i = 0; result && (i < DIM); i++)
      result = result && (d_vector[i] > rhs.d_vector[i]);

   return result;
}

template<int DIM> inline
bool IntVector<DIM>::operator>=(const IntVector<DIM>& rhs) const
{
   bool result = true;
   
   for(int i = 0; result && (i < DIM); i++)
      result = result && (d_vector[i] >= rhs.d_vector[i]);

   return result;
}

template<int DIM> inline
int IntVector<DIM>::min() const
{
   int min = d_vector[0];

   for(int i = 1; i < DIM; i++)
      if (d_vector[i] < min) 
	 min = d_vector[i];

   return(min);
}

template<int DIM> inline
void IntVector<DIM>::max(const IntVector<DIM>& rhs)
{
   for(int i = 0; i < DIM; i++)
      if (rhs.d_vector[i] > d_vector[i]) 
	 d_vector[i] = rhs.d_vector[i];
}

template<int DIM> inline
int IntVector<DIM>::max() const
{
   int max = d_vector[0];

   for(int i = 1; i < DIM; i++)
      if (d_vector[i] > max) 
	 max = d_vector[i];

   return(max);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::min(const IntVector<DIM>& a,
                                     const IntVector<DIM>& b)
{
   IntVector<DIM> tmp = a;
   tmp.min(b);
   return(tmp);
}

template<int DIM> inline
IntVector<DIM> IntVector<DIM>::max(const IntVector<DIM>& a,
                                     const IntVector<DIM>& b)
{
   IntVector<DIM> tmp = a;
   tmp.max(b);
   return(tmp);
}

template<int DIM> inline
int IntVector<DIM>::getProduct() const
{
   int prod = 1;
   for (int i = 0; i < DIM; i++) {
     prod *= d_vector[i];
   }
   return(prod);
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
