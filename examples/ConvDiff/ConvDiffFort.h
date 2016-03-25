//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/ConvDiff/ConvDiffFort.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
