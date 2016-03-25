/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/printObject.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Misc printing functions in FAC solver test.
 */

#ifndef included_printObject_h
#define included_printObject_h


#include "tbox/Pointer.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchData.h"
#include "ArrayData.h"
#include "Box.h"

#include <string>

using namespace std;
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
