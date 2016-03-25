/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/QuarticFcn.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description:	Quartic function functor.
 */

#ifndef included_QuarticFcn_h
#define included_QuarticFcn_h

#include <iostream>

using namespace std;


/*!
  @brief Quartic function functor.

  This functor exists for convenience, not efficiency.
  If you use it heavily, you can improve your efficiency
  by replacing its usage with lower level codes.
*/
class QuarticFcn {

public:

  /*!
    @brief Default constructor

    Construct a trivial polynomial function.
  */
  QuarticFcn();

  /*!
    @brief Copy constructor
  */
  QuarticFcn( const QuarticFcn &r );

  /*!
    @brief Mapping from polynomial term to coefficient number.

    Each member of this enumeration corresponds to a term in
    the polynomial.  @c co_c corresponds to the constant term.
    The rest should be self-explanatory.
  */
  enum ctype {
    co_c	= 0
#if NDIM >= 1
    , co_x	= 1
    , co_xx	= 2
    , co_xxx	= 3
    , co_xxxx	= 4
#endif
#if NDIM >= 2
    , co_y	= 5
    , co_xy	= 6
    , co_yy	= 7
    , co_xxy	= 8
    , co_xyy	= 9
    , co_yyy	= 10
    , co_xxxy	= 11
    , co_xxyy	= 12
    , co_xyyy	= 13
    , co_yyyy	= 14
#endif
#if NDIM >= 3
    , co_z	= 15
    , co_xz	= 16
    , co_yz	= 17
    , co_zz	= 18
    , co_xxz	= 19
    , co_xyz	= 20
    , co_xzz	= 21
    , co_yyz	= 22
    , co_yzz	= 23
    , co_zzz	= 24
    , co_xxxz	= 25
    , co_xxyz	= 26
    , co_xxzz	= 27
    , co_xyyz	= 28
    , co_xyzz	= 29
    , co_xzzz	= 30
    , co_yyyz	= 31
    , co_yyzz	= 32
    , co_yzzz	= 33
    , co_zzzz	= 34
#endif
  };

  /*!
    @brief Set polynomial coefficient corresponding to indicated term.
  */
  int setPolynomialCoef( /*! Term describer */ ctype type,
			 /*! Coefficient */ double c );

  /*!
    @brief Access to polynomial coefficient.
  */
  double getPolynomialCoef( /*! Term describer */ ctype type ) const;

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


  //@{

  //! @name Arithmetic operators.

  /*!
    @brief Addition of two polynomials.
  */
  QuarticFcn operator+( const QuarticFcn &r ) const;

  /*!
    @brief Subtraction of two polynomials.
  */
  QuarticFcn operator-( const QuarticFcn &r ) const;

  /*!
    @brief Equal operator.
  */
  QuarticFcn &operator=( const QuarticFcn &r );

  //@}


  //@{

  //! @name Differential operations.

  /*!
    @brief Differentiate and return new function.

    @return Differentiated function.
  */
QuarticFcn differentiate( unsigned short int x
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
QuarticFcn &differentiateSelf( unsigned short int x
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

    @see QuarticFcn::ctype
  */
  friend std::istream &operator>>( std::istream &ci, QuarticFcn &cf );
  /*!
    @brief Output to stream.

    Outputs sets of integer=double values where the integer indicates
    the coefficient describer (see ctype) and the double is the
    value of the coefficient.

    @see QuarticFcn::ctype
  */
  friend std::ostream &operator<<( std::ostream &co, const QuarticFcn &cf );
  //@}

  /*!
    @brief Dimensionally-dependent number of polynomial coefficients.

    The number of coefficients are the number of ctype enumerations.
  */
#if NDIM == 1
  #define NUMBER_OF_COEF 5
#endif
#if NDIM == 2
  #define NUMBER_OF_COEF 15
#endif
#if NDIM == 3
  #define NUMBER_OF_COEF 35
#endif


private:

  /*!
    @brief tbox::List of coefficients.

    The ordering is specified in the enum.
  */
  double d_coefs[NUMBER_OF_COEF];

};




#endif
