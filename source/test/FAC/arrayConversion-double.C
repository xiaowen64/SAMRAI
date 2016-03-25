/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/arrayConversion-double.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

