/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/setArrayData.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Misc array setting functions in FAC solver test.
 */

#include "SAMRAI_config.h"

#include "setArrayData.h"
#include "MDA_Access.h"
#include <math.h>

#include <iostream>

using namespace std;

/*!
  \file
  \brief Set array data when array is given as a pointer or
  MultiDimArrayAccess object.

  There is no SAMRAI-specific code in this file.
*/

void setArrayDataToConstant(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower, const int *upper
  , const double *xlo, const double *xhi, const double *h
  , double value
) {
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      s(i,j) = value;
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	s(i,j,k) = value;
      }
    }
  }
#endif
  return;
}

void setArrayDataTo(
  double *ptr
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *coef
) {
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > s( ptr , lower , upper );
  setArrayDataTo( s, lower, upper, xlo, xhi, h, coef );
  return;
}

void setArrayDataTo(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *coef
) {
#if NDIM == 2
  const double ucoef[NDIM] = { 1., 1.};
  if ( coef == NULL ) coef = ucoef;
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      s(i,j) = coef[0]*x + coef[1]*y;
    }
  }
#endif
#if NDIM == 3
  const double ucoef[NDIM] = { 1., 1., 1. };
  if ( coef == NULL ) coef = ucoef;
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	s(i,j,k) = coef[0]*x + coef[1]*y + coef[2]*z;
      }
    }
  }
#endif
  return;
}

void setArrayDataToSinusoidal(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
  , const double *npi, const double *ppi
) {
  double nx=npi[0], px=ppi[0];
#if NDIM >= 2
  double ny=npi[1], py=ppi[1];
#endif
#if NDIM >= 3
  double nz=npi[2], pz=ppi[2];
#endif
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    double siny = sin(M_PI*(ny*y+py));
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      double sinx = sin(M_PI*(nx*x+px));
      s(i,j) = sinx * siny;
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    double sinz = sin(M_PI*(nz*z+pz));
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      double siny = sin(M_PI*(ny*y+py));
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	double sinx = sin(M_PI*(nx*x+px));
	s(i,j,k) = sinx * siny * sinz;
      }
    }
  }
#endif
  return;
}

void setArrayDataToSinusoidalGradient(
    double **g_ptr
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
) {
  double *gx_ptr = g_ptr[0];
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > gx( gx_ptr , lower , upper );
#if NDIM >= 2
  double *gy_ptr = g_ptr[1];
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > gy( gy_ptr , lower , upper );
#endif
#if NDIM >= 3
  double *gz_ptr = g_ptr[2];
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > gz( gz_ptr , lower , upper );
#endif
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    double siny = sin(2*M_PI*y);
    double cosy = cos(2*M_PI*y);
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      double sinx = sin(2*M_PI*x);
      double cosx = cos(2*M_PI*x);
      gx(i,j) = 2*M_PI*cosx * siny;
      gy(i,j) = sinx * 2*M_PI*cosy;
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    double sinz = sin(2*M_PI*z);
    double cosz = cos(2*M_PI*z);
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      double siny = sin(2*M_PI*y);
      double cosy = cos(2*M_PI*y);
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	double sinx = sin(2*M_PI*x);
	double cosx = cos(2*M_PI*x);
	gx(i,j,k) = 2*M_PI*cosx * siny * sinz;
	gy(i,j,k) = sinx * 2*M_PI*cosy * sinz;
	gz(i,j,k) = sinx * cosy * 2*M_PI*cosz;
      }
    }
  }
#endif
}

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
) {
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      s(i,j) = a0 + ax*x + ay*y + axy*x*y;
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	s(i,j,k) = a0 + ax*x + ay*y + az*z
	  + axy*x*y + axz*x*z + ayz*y*z + axyz*x*y*z;
      }
    }
  }
#endif
  return;
}

void setArrayDataToScaled(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower, const int *upper
  , double factor
) {
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      s(i,j) *= factor;
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	s(i,j,k) *= factor;
      }
    }
  }
#endif
  return;
}

/*!
  \brief Set array to Michael's exact solution
*/
void setArrayDataToPerniceExact(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
) {
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      s(i,j) = x*(1-x)*y*(1-y);
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	s(i,j,k) = x*(1-x)*y*(1-y)*z*(1-z);
      }
    }
  }
#endif
  return;
}

/*!
  \brief Set array to Michael's source
*/
void setArrayDataToPerniceSource(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *lower
  , const int *upper
  , const double *xlo, const double *xhi, const double *h
) {
#if NDIM == 2
  for ( int j=lower[1]; j<=upper[1]; ++j) {
    double y = xlo[1] + h[1]*(j-lower[1]+0.5);
    for ( int i=lower[0]; i<=upper[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-lower[0]+0.5);
      s(i,j) = -2*( x*(1-x) + y*(1-y) );
    }
  }
#endif
#if NDIM == 3
  for ( int k=lower[2]; k<=upper[2]; ++k) {
    double z = xlo[2] + h[2]*(k-lower[2]+0.5);
    for ( int j=lower[1]; j<=upper[1]; ++j) {
      double y = xlo[1] + h[1]*(j-lower[1]+0.5);
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-lower[0]+0.5);
	s(i,j,k) = -2*( x*(1-x) + y*(1-y) + z*(1-z) );
      }
    }
  }
#endif
  return;
}

void setArrayDataToSinusoid(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *beg
  , const int *end
  , const int *ilo, const double *xlo, const double *h
  , const SinusoidFcn &fcn )
{
#if NDIM == 2
  for ( int j=beg[1]; j<=end[1]; ++j) {
    double y = xlo[1] + h[1]*(j-ilo[1]+0.5);
    for ( int i=beg[0]; i<=end[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-ilo[0]+0.5);
      s(i,j) = fcn(x,y);
    }
  }
#endif
#if NDIM == 3
  for ( int k=beg[2]; k<=end[2]; ++k) {
    double z = xlo[2] + h[2]*(k-ilo[2]+0.5);
    for ( int j=beg[1]; j<=end[1]; ++j) {
      double y = xlo[1] + h[1]*(j-ilo[1]+0.5);
      for ( int i=beg[0]; i<=end[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-ilo[0]+0.5);
	s(i,j,k) = fcn(x,y,z);
      }
    }
  }
#endif
  return;
}

void setArrayDataToQuartic(
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > &s
  , const int *beg
  , const int *end
  , const int *ilo, const double *xlo, const double *h
  , const QuarticFcn &fcn )
{
#if NDIM == 2
  for ( int j=beg[1]; j<=end[1]; ++j) {
    double y = xlo[1] + h[1]*(j-ilo[1]+0.5);
    for ( int i=beg[0]; i<=end[0]; ++i ) {
      double x = xlo[0] + h[0]*(i-ilo[0]+0.5);
      s(i,j) = fcn(x,y);
    }
  }
#endif
#if NDIM == 3
  for ( int k=beg[2]; k<=end[2]; ++k) {
    double z = xlo[2] + h[2]*(k-beg[2]+0.5);
    for ( int j=beg[1]; j<=end[1]; ++j) {
      double y = xlo[1] + h[1]*(j-ilo[1]+0.5);
      for ( int i=beg[0]; i<=end[0]; ++i ) {
	double x = xlo[0] + h[0]*(i-ilo[0]+0.5);
	s(i,j,k) = fcn(x,y,z);
      }
    }
  }
#endif
  return;
}
