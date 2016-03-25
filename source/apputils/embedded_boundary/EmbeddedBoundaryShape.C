//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/embedded_boundary/EmbeddedBoundaryShape.C $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
   NULL_USE(xyz);
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
   NULL_USE(nx);
   NULL_USE(dx);
   NULL_USE(origin);
   NULL_USE(inout);

   TBOX_ERROR("EmbeddedBoundaryShape::isInside(): "
                << "\nNo implementation provided for this shape. "
                << std::endl);
}
        
}
}
#endif
