/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/FAC/GaussianFcn.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2172 $
 * Modified:    $LastChangedDate: 2008-05-02 11:02:08 -0700 (Fri, 02 May 2008) $
 * Description: Gaussian function support for FAC solver tests.
 */

#include "SAMRAI_config.h"

#include "GaussianFcn.h"
#include <math.h>
#include <stdlib.h>
#include "tbox/Utilities.h"

#include <string>
#include <string.h>

/*
  Temporary fix for g++ lacking instantiations when --no-implicit-templates
  is used (by SAMRAI)
*/
#define fill_n(p,n,v) { size_t _i; for ( _i=0; _i<n; ++_i ) p[_i]=v; }
#define copy_n(s,n,d) { size_t _i; for ( _i=0; _i<n; ++_i ) d[_i]=s[_i]; }

GaussianFcn::GaussianFcn() 
: d_amp(1.0) ,
  d_lambda(-1.0)
{
  fill_n( d_center, NDIM, 0.0 )
}

int GaussianFcn::setAmplitude( const double amp ) {
  d_amp = amp;
  return 0;
}

int GaussianFcn::setLambda( const double lambda ) {
  d_lambda = lambda;
  return 0;
}

int GaussianFcn::setCenter( const double *center ) {
  for ( size_t i=0; i<NDIM; ++i ) d_center[i] = center[i];
  return 0;
}

double GaussianFcn::getAmplitude() const {
  return d_amp;
}

double GaussianFcn::getLambda() const {
  return d_lambda;
}

int GaussianFcn::getCenter( double *center ) const {
  for ( size_t i=0; i<NDIM; ++i ) center[i] = d_center[i];
  return 0;
}

#if NDIM == 1
double GaussianFcn::operator()( double x ) const {
  double rval;
  rval = (x-d_center[0])*(x-d_center[0]);
  rval = exp( d_lambda*rval );
  return rval;
}
#endif
#if NDIM == 2
double GaussianFcn::operator()( double x, double y ) const {
  double rval;
  rval = (x-d_center[0])*(x-d_center[0]) + (y-d_center[1])*(y-d_center[1]);
  rval = exp( d_lambda*rval );
  return rval;
}
#endif
#if NDIM == 3
double GaussianFcn::operator()( double x, double y, double z ) const {
  double rval;
  rval =
    (x-d_center[0])*(x-d_center[0]) +
    (y-d_center[1])*(y-d_center[1]) +
    (z-d_center[2])*(z-d_center[2]) ;
  rval = exp( d_lambda*rval );
  return rval;
}
#endif

#define EAT_WS(s) { while ( s.peek() == ' '			\
			 || s.peek() == '\t'			\
                         || s.peek() == '\n' ) { s.get(); } }

istream &operator>>( istream &ci, GaussianFcn &gf ) {
  fill_n( gf.d_center, NDIM, 0.0 )
  gf.d_amp = 1.0;
  gf.d_lambda = -1.0;
  char dummy, name[6];
  EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
  ci >> dummy;
  TBOX_ASSERT( dummy == '{' );
  EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
  while ( ci.peek() != '}' ) {
    ci.read(name,2);
    if ( name[0]=='l' ) {
      // Expect form lambda=<float>
      ci.read(name,5);
      TBOX_ASSERT( !strncmp(name,"mbda=",5) );
      ci >> gf.d_lambda;
      EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
    }
    else if ( name[0]=='a' ) {
      // Expect form amp=<float>
      ci.read(name,2);
      TBOX_ASSERT( !strncmp(name,"p=",2) );
      ci >> gf.d_amp;
      EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
    }
    else if ( name[0]=='c' ) {
      // Expect form c[xyz]=<float>
      unsigned short dim( name[1]=='x' ? 0 :
			  name[1]=='y' ? 1 :
			  name[1]=='z' ? 2 : 3 );
      TBOX_ASSERT( dim < NDIM );
      ci >> dummy; TBOX_ASSERT( dummy == '=' );
      ci >> gf.d_center[dim];
      EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
    }
    else {
      abort();
    }
  }
  return ci;
}

ostream &operator<<( ostream &co, const GaussianFcn &gf ) {
  co << "{ amp=" << gf.d_amp << " lambda=" << gf.d_lambda
     << " cx=" << gf.d_center[0];
#if NDIM >= 2
  co << " cy=" << gf.d_center[1];
#endif
#if NDIM >= 3
  co << " cz=" << gf.d_center[2];
#endif
  co << " }";
  return co;
}
