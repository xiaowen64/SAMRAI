#ifndef included_SinusoidFcn_h
#define included_SinusoidFcn_h

#include <iostream>

using namespace std;

/*!
  @brief Sinusoid function functor.

  This functor exists for convenience, not efficiency.
  If you use it heavily, you can improve your efficiency
  by replacing its usage with lower level codes.
*/
class SinusoidFcn {

public:

  SinusoidFcn();


  /*!
    @brief Return the function value.
  */
#if NDIM == 1
  double operator()( double x ) const;
#endif
#if NDIM == 2
  double operator()( double x, double y ) const;
#endif
#if NDIM == 3
  double operator()( double x, double y, double z ) const;
#endif


  /*!
    @brief Set amplitude.
  */
  int setAmplitude( const double amp );

  /*!
    @brief Set all wave numbers.

    Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int setWaveNumbers( const double *npi );

  /*!
    @brief Set all phase angles.

    Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int setPhaseAngles( const double *ppi );


  /*!
    @brief Get wave numbers.

    Wave numbers are be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int getWaveNumbers( double *npi ) const;

  /*!
    @brief Get phase angles.

    Wave numbers are be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int getPhaseAngles( double *ppi ) const;


  //@{

  //! @name Differential operations.

  /*!
    @brief Differentiate and return new function.

    @return Differentiated function.
  */
SinusoidFcn differentiate( unsigned short int x
#if NDIM >= 2
			    , unsigned short int y=0
#endif
#if NDIM >= 3
			    , unsigned short int z=0
#endif
			    ) const;

  /*!
    @brief Differentiate self and return self reference.

    @return Self reference
  */
SinusoidFcn &differentiateSelf( unsigned short int x
#if NDIM >= 2
				 , unsigned short int y=0
#endif
#if NDIM >= 3
				 , unsigned short int z=0
#endif
				 );

  //@}


  //@{
  /*!
    @name IO operators.
  */
  /*!
    @brief Input from stream.

    Stream extract function reads input in the format used
    by the stream insert function.  Except for allowing for
    missing coefficients (set to zero), this function requires
    the correct syntax or the result is undefined.

    @see SinusoidFcn::ctype
  */
  friend istream &operator>>( istream &ci, SinusoidFcn &cf );
  /*!
    @brief Output to stream.

    Outputs sets of name=double values where the name is one
    of nx, px, ny, nz or pz, and the double is the value of
    the coefficient.

    @see SinusoidFcn::ctype
  */
  friend ostream &operator<<( ostream &co, const SinusoidFcn &cf );
  //@}


private:

  //! Amplitude
  double d_amp;
  //! Wave number in half-cycles
  double d_npi[NDIM];
  //! Phase shift in half-cycles
  double d_ppi[NDIM];

};




#endif
