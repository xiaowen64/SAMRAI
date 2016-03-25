//
// File:        EmbeddedBoundaryShape.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 605 $
// Modified:    $Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description: Base class for analytic embedded Boundaries
//              
// 

#ifndef included_appu_EmbeddedBoundaryShape_C
#define included_appu_EmbeddedBoundaryShape_C

#include "EmbeddedBoundaryShape.h"

#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
   namespace appu {


/*
*******************************************************************
* 
*  Empty constructor and destructor
*
*******************************************************************
*/
template<int DIM> EmbeddedBoundaryShape<DIM>::EmbeddedBoundaryShape()
{
}

template<int DIM> EmbeddedBoundaryShape<DIM>::~EmbeddedBoundaryShape<DIM>()
{   
}

/*
*************************************************************************
*                                                                       *
* Default virtual function implementations.                             *
*                                                                       *
*************************************************************************
*/
template<int DIM> bool 
EmbeddedBoundaryShape<DIM>::isInside(const double* xyz) const
{
   (void) xyz;
   TBOX_ERROR("EmbeddedBoundaryShape::isInside(): "
                << "\nNo implementation provided for this shape."
                << endl);
   return(false);
}


template<int DIM> void 
EmbeddedBoundaryShape<DIM>::isInside(const int* nx,
                                     const double* dx,
                                     const double* origin,
                                     int* inout) const
{
   (void) nx;
   (void) dx;
   (void) origin;
   (void) inout;

   TBOX_ERROR("EmbeddedBoundaryShape::isInside(): "
                << "\nNo implementation provided for this shape. "
                << endl);
}
        
}
}
#endif
