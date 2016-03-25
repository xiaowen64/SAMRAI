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

