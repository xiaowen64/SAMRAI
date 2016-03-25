/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/FAC/printObject.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Misc printing functions in FAC solver test.
 */

#include "SAMRAI_config.h"

#include "printObject.h"

#include "MDA_Access.h"

#include "SAMRAI_config.h"
#include "tbox/Pointer.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchData.h"
#include "Box.h"
#include "ArrayData.h"

using namespace SAMRAI;

/*!
  \brief Print a box
*/
int printObject(
  ostream &os ,
  const hier::Box<NDIM> &box ,
  const string& border ,
  unsigned short depth )
{
#if NDIM == 1
  os << border << "( " << box.lower(0)
     << " ) to ( "
     << box.upper(0)
     << " )\n";
#endif
#if NDIM == 2
  os << border << "( " << box.lower(0) << ' ' << box.lower(1)
     << " ) to ( "
     << box.upper(0) << ' ' << box.upper(1)
     << " )\n";
#endif
#if NDIM == 3
  os << border << "( "
     << box.lower(0) << ' ' << box.lower(1) << ' ' << box.lower(2)
     << " ) to ( "
     << box.upper(0) << ' ' << box.upper(1) << ' ' << box.upper(2)
     << " )\n";
#endif
  return 0;
}

/*!
  \brief Print a patch data object
*/
int printObject(
  ostream &os ,
  const hier::PatchData<NDIM> &pdat ,
  const string& border ,
  unsigned short depth )
{
  const hier::Box<NDIM> &rbox = pdat.getBox();
  const hier::Box<NDIM> &gbox = pdat.getBox();
  os << border
     << "Box = ( "
     << rbox.lower(0) << ' ' << rbox.lower(1)
     << " ) to ( "
     << rbox.upper(0) << ' ' << rbox.upper(1)
     << " )\n";
  os << border
     << "Ghost hier::Box = ( "
     << gbox.lower(0) << ' ' << gbox.lower(1)
     << " ) to ( "
     << gbox.upper(0) << ' ' << gbox.upper(1)
     << " )\n";
  return 0;
}

/*!
  \brief Print an array data object
*/
template <class T>
int printObject(
  ostream &os ,
  const pdat::ArrayData<NDIM,T> &adat ,
  const int depth ,
  const string& border )
{
  os << border << "printObject(pdat::ArayData)... " << endl;
  const hier::Box<NDIM> &rbox = adat.getBox();
  os << border
     << "Box = ( "
     << rbox.lower(0) << ' ' << rbox.lower(1)
     << " ) to ( "
     << rbox.upper(0) << ' ' << rbox.upper(1)
     << " )\n";
  os << border << "depth " << depth << "\n";
  printObject( os
	     , (const T*)adat.getPointer(depth)
	     , (const int*)adat.getBox().lower()
	     , (const int*)adat.getBox().upper()
	     );
  return 0;
}
template int printObject<double>( ostream &os , const pdat::ArrayData<NDIM,double> &adat, const int depth, const string& border );

/*!
  \brief Print an array
*/
int printObject(
  ostream &os
, const double *a_ptr , const int *a_lower , const int *a_upper
, const int *lower
, const int *upper
) {
  MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > a( a_ptr , a_lower , a_upper );
  if ( ! lower ) lower = a_lower;
  if ( ! upper ) upper = a_upper;
  os << "\narray data... " << endl;
  for ( int j=lower[1]; j<=upper[1]; ++j ) {
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      os << i << ' ' << j << ' ' << a(i,j) << endl;
    }
  }
  return 0;
}

