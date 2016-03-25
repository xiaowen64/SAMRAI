


#ifndef included_printObject_h
#define included_printObject_h


#include "tbox/Pointer.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchData.h"
#include "ArrayData.h"
#include "Box.h"

#include <string>

using namespace SAMRAI;


/*!
  @brief Print a box
*/
int printObject(
  ostream &os ,
  const hier::Box<NDIM> &box ,
  const string& border=string() ,
  unsigned short depth=0 );


/*!
  @brief Print a patch data object
*/
int printObject(
  ostream &os ,
  const hier::PatchData<NDIM> &pdat ,
  const string& border=string() ,
  unsigned short depth=0 );


/*!
  @brief Print an array
*/
int printObject(
  ostream &os ,
  const double *a_ptr , const int *a_lower , const int *a_upper ,
  const int *lower=NULL ,
  const int *upper=NULL );



/*!
  @brief Print an array data object
*/
template <class T>
int printObject(
  ostream &os ,
  const pdat::ArrayData<NDIM,T> &adat ,
  const int depth=0 ,
  const string& border=string() );


#endif	// included_printObject_h
