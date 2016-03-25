//
// File:        HeatEqnFort.h
// Package:     SAMRAI application
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: F77 external declarations for SAMRAI Heat Equation example.
//

#include <math.h>
#include <signal.h>

// Link between C/C++ and Fortran files
//       name in             name in
//      C/C++ code            Fortran code
//      ----------            ------------
#define FORT_SPHERE_INIT      initsphere_
#define FORT_COMP_RHS         computerhs_
#define FORT_RK_STEP          rkstep_
#define FORT_TAG_CELLS        tagcells_



// Function argument list interfaces
extern "C" {
  void FORT_SPHERE_INIT(
  const double*, const double*, const double*, 
  const int& , const int& , const int& , const int& , 
#if (NDIM==3)
  const int& , const int& , 
#endif
  const int& , const int& , 
#if (NDIM==3)
  const int& , 
#endif
  double*,
  const double*,
  const double*,
  const double*,
  const double&,
  const int&);

  void FORT_COMP_RHS(
  const int& , const int& , const int& , const int& , 
#if (NDIM==3)
  const int& , const int& , 
#endif
  const int& , const int& , 
#if (NDIM==3)
  const int& , 
#endif
  const double*,  // dx
  const double*,  // d_convection_coeff
  const double&,  // d_diffusion_coeff
  const double&,  // d_source_coeff
  double*,        // prim_var_updated
  double*,        // function_eval
  const int&);    // NEQU

  void FORT_RK_STEP(
  const int& , const int& , const int& , const int& , 
#if (NDIM==3)
  const int& , const int& , 
#endif
  const int& , const int& , 
#if (NDIM==3)
  const int& , 
#endif
  const double&, const double&, const double&, const double&,
  const double*, 
  const double&, 
  const double&,
  double*,
  const double*,
  const double*,
  const int&);


  void FORT_TAG_CELLS( 
  const int& , const int& , const int& , const int& , 
#if (NDIM==3)
  const int& , const int& , 
#endif
  int*, 
  const double*,
  const int&,
  const double*,
  const int&);


  void FORT_DIRICHLET_BCS( const int&,
  const int& , const int& , const int& , const int& ,
#if (NDIM==3)
  const int& , const int& ,
#endif
  const int& , const int& , const int& , const int& ,
#if (NDIM==3)
  const int& , const int& ,
#endif
  const int& , const int& ,  // d_nghosts(0), d_nghosts(1)
#if (NDIM==3)
  const int& ,               // d_nghosts(2)
#endif
  const int& ,               // bdry_index
  double*,                   // primitive vars
  const double*,             // bdry_values
  const int&);               // NEQU

  void FORT_NEUMANN_BCS( const int&,
  const int& , const int& , const int& , const int& ,
#if (NDIM==3)
  const int& , const int& ,
#endif
  const int& , const int& , const int& , const int& ,
#if (NDIM==3)
  const int& , const int& ,
#endif
  const int& , const int& ,  // d_nghosts(0), d_nghosts(1)
#if (NDIM==3)
  const int& ,               // d_nghosts(2)
#endif
  const int& ,               // bdry_index
  double*,                   // primitive vars
  const double*,             // dx
  const double*,             // bdry_values
  const int&);               // NEQU

}
