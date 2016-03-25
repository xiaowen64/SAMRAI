/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	Sinusoid function functor.
*/

#include "SAMRAI_config.h"

#include "SinusoidFcn.h"
#include <math.h>
#include <assert.h>

/*
  Temporary fix for g++ lacking instantiations when --no-implicit-templates
  is used (by SAMRAI)
*/
#define fill_n(p,n,v) { size_t _i; for ( _i=0; _i<n; ++_i ) p[_i]=v; }
#define copy_n(s,n,d) { size_t _i; for ( _i=0; _i<n; ++_i ) d[_i]=s[_i]; }

SinusoidFcn::SinusoidFcn() 
: d_amp(1.0)
{
  fill_n( d_npi, NDIM, 0.0 )
  fill_n( d_ppi, NDIM, 0.0 )
}

int SinusoidFcn::setAmplitude( const double amp ) {
  d_amp = amp;
  return 0;
}

int SinusoidFcn::setWaveNumbers( const double *npi ) {
  for ( size_t i=0; i<NDIM; ++i ) d_npi[i] = npi[i];
  return 0;
}

int SinusoidFcn::getWaveNumbers( double *npi ) const {
  for ( size_t i=0; i<NDIM; ++i ) npi[i] = d_npi[i];
  return 0;
}

int SinusoidFcn::setPhaseAngles( const double *ppi ) {
  for ( size_t i=0; i<NDIM; ++i ) d_ppi[i] = ppi[i];
  return 0;
}

int SinusoidFcn::getPhaseAngles( double *ppi ) const {
  for ( size_t i=0; i<NDIM; ++i ) ppi[i] = d_ppi[i];
  return 0;
}

#if NDIM == 1
double SinusoidFcn::operator()( double x ) const {
  double rval;
  rval = d_amp
       * sin(M_PI*(d_npi[0]*x+d_ppi[0]));
  return rval;
}
#endif
#if NDIM == 2
double SinusoidFcn::operator()( double x, double y ) const {
  double rval;
  rval = d_amp
       * sin(M_PI*(d_npi[0]*x+d_ppi[0]))
       * sin(M_PI*(d_npi[1]*y+d_ppi[1]));
  return rval;
}
#endif
#if NDIM == 3
double SinusoidFcn::operator()( double x, double y, double z ) const {
  double rval;
  rval = d_amp
       * sin(M_PI*(d_npi[0]*x+d_ppi[0]))
       * sin(M_PI*(d_npi[1]*y+d_ppi[1]))
       * sin(M_PI*(d_npi[2]*z+d_ppi[2]));
  return rval;
}
#endif

SinusoidFcn &SinusoidFcn::differentiateSelf( unsigned short int x
#if NDIM >= 2
					     , unsigned short int y
#endif
#if NDIM >= 3
					     , unsigned short int z
#endif
					     )
{
  /*
    Since differentiation commutes,
    simply differentiate one direction at a time.
  */
#if NDIM >= 1
  // Differentiate in x direction.
  for ( ; x>0; --x ) {
    d_amp *= M_PI*d_npi[0];
    d_ppi[0] += 0.5;
  }
  while ( d_ppi[0] > 2 ) d_ppi[0] -= 2.0;
#endif
#if NDIM >= 2
  // Differentiate in y direction.
  for ( ; y>0; --y ) {
    d_amp *= M_PI*d_npi[1];
    d_ppi[1] += 0.5;
  }
  while ( d_ppi[1] > 2 ) d_ppi[1] -= 2.0;
#endif
#if NDIM >= 3
  // Differentiate in z direction.
  for ( ; z>0; --z ) {
    d_amp *= M_PI*d_npi[2];
    d_ppi[2] += 0.5;
  }
  while ( d_ppi[2] > 2 ) d_ppi[2] -= 2.0;
#endif
  return *this;
}

SinusoidFcn SinusoidFcn::differentiate( unsigned short int x
#if NDIM >= 2
					, unsigned short int y
#endif
#if NDIM >= 3
					, unsigned short int z
#endif
					) const
{
  SinusoidFcn rval(*this);
  rval.differentiateSelf( x
#if NDIM >= 2
			, y
#endif
#if NDIM >= 3
			, z
#endif
			);
  return rval;
}

#define EAT_WS(s) { while ( s.peek() == ' '			\
			 || s.peek() == '\t'			\
                         || s.peek() == '\n' ) { s.get(); } }

istream &operator>>( istream &ci, SinusoidFcn &sf ) {
  fill_n( sf.d_npi, NDIM, 0.0 )
  fill_n( sf.d_ppi, NDIM, 0.0 )
  char dummy, name[2];
  EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
  ci >> dummy;
  assert ( dummy == '{' );
  EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
  while ( ci.peek() != '}' ) {
    ci.read(name,2);
    if ( name[0] == 'a' && name[1] == 'm' ) {
      ci >> dummy; assert ( dummy == 'p' );
      ci >> dummy; assert ( dummy == '=' );
      ci >> sf.d_amp;
    }
    else {
      ci >> dummy; assert( dummy == '=' );
      double *data( name[0]=='n' ? sf.d_npi : sf.d_ppi );
      unsigned short dim( name[1]=='x' ? 0 :
			  name[1]=='y' ? 1 :
			  name[1]=='z' ? 2 : 3 );
      assert( dim < NDIM );
      ci >> data[dim];
    }
    EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
  }
  return ci;
}

ostream &operator<<( ostream &co, const SinusoidFcn &sf ) {
  co << "{ amp=" << sf.d_amp;
  co << " nx=" << sf.d_npi[0] << " px=" << sf.d_npi[0];
#if NDIM >= 2
  co << " ny=" << sf.d_npi[1] << " py=" << sf.d_npi[1];
#endif
#if NDIM >= 3
  co << " nz=" << sf.d_npi[2] << " pz=" << sf.d_npi[2];
#endif
  co << " }";
  return co;
}
