/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/arrayConversion.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Array conversion functions in FAC solver test.
 */

#ifndef included_arrayConversion_h
#define included_arrayConversion_h

#include "MDA_Access.h"
#include "ArrayData.h"

using namespace SAMRAI;

/*!
  @brief Generate a MultiDimArrayAccess object
  from a pdat::ArrayData<NDIM> object.
*/
template < class T >
MDA_Access<T,NDIM,MDA_OrderColMajor<NDIM> > arrayData2ArrayAccess (
  pdat::ArrayData<NDIM,T> &adat ,
  int depth=0 )
{
  return MDA_Access<T,NDIM,MDA_OrderColMajor<NDIM> >(
	   adat.getPointer(depth) ,
	   (const int*)adat.getBox().lower() ,
	   (const int*)adat.getBox().upper() );
}

/*!
  @brief Generate a const MultiDimArrayAccess object
  from a const pdat::ArrayData<NDIM> object.
*/
template < class T >
MDA_AccessConst<T,NDIM,MDA_OrderColMajor<NDIM> > arrayData2ArrayAccess (
  const pdat::ArrayData<NDIM,T> &adat ,
  int depth=0 )
{
  return MDA_AccessConst<T,NDIM,MDA_OrderColMajor<NDIM> >(
	   adat.getPointer(depth) ,
	   (const int*)adat.getBox().lower() ,
	   (const int*)adat.getBox().upper() );
}

#endif	// include_arrayConversion_h
