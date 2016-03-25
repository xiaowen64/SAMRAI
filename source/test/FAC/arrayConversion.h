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
