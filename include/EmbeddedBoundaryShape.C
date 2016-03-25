//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/embedded_boundary/EmbeddedBoundaryShape.C $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Base class for analytic embedded Boundaries
//              
// 

#ifndef included_appu_EmbeddedBoundaryShape_C
#define included_appu_EmbeddedBoundaryShape_C

#include "EmbeddedBoundaryShape.h"

#include "tbox/Utilities.h"


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
                << std::endl);
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
                << std::endl);
}
        
}
}
#endif
