/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/setArrayData.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Misc array setting functions in FAC solver test.
 */

#ifndef include_setArrayData_h
#define include_setArrayData_h

#include "MDA_Access.h"
#include "QuarticFcn.h"
#include "SinusoidFcn.h"



void setArrayDataTo(
  double *ptr
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *coef=NULL
);

void setArrayDataTo(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *coef=NULL
);

void setArrayDataToSinusoidal(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *npi, const double *ppi
);

void setArrayDataToSinusoidalGradient(
    double **g_ptr
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
);

void setArrayDataToConstant(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
, const int *lower, const int *upper
, const double *xlo , const double *xhi , const double *h
, double value
);

void setArrayDataToLinear(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s,
  const int *lower,
  const int *upper,
  const double *xlo, const double *xhi, const double *h,
#if NDIM == 2
  double a0, double ax, double ay, double axy
#endif
#if NDIM == 3
  double a0, double ax, double ay, double az,
  double axy, double axz, double ayz, double axyz
#endif
);

void setArrayDataToScaled(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower, const int *upper
  , double factor
);

void setArrayDataToPerniceExact(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
, const int *lower
, const int *upper
, const double *xlo , const double *xhi , const double *h
);

void setArrayDataToPerniceSource(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
, const int *lower
, const int *upper
, const double *xlo , const double *xhi , const double *h
);
void setArrayDataToSinusoid(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *beg
  , const int *end
  , const int *ilo, const double *xlo, const double *h
  , const SinusoidFcn &fcn );
void setArrayDataToQuartic(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *beg
  , const int *end
  , const int *ilo, const double *xlo, const double *h
  , const QuarticFcn &fcn );

#endif	// include_setArrayData_h
