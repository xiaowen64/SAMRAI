/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/FAC/arrayConversion-double.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Instantiate array conversion functions in FAC solver test.
 */

#include "SAMRAI_config.h"
#include "arrayConversion.h"

template
MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > arrayData2ArrayAccess (
  pdat::ArrayData<NDIM,double> &adat ,
  int depth );

template
MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > arrayData2ArrayAccess (
  const pdat::ArrayData<NDIM,double> &adat ,
  int depth );

