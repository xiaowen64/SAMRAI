//
// File:        EmbeddedBoundaryShapeSphere.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 605 $
// Modified:    $Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description: Sphere embedded boundary shape
//              
// 

#ifndef included_appu_EmbeddedBoundaryShapeSphere_C
#define included_appu_EmbeddedBoundaryShapeSphere_C

#include "EmbeddedBoundaryShapeSphere.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#define included_assert
#endif

#include "tbox/IEEE.h"

#ifdef DEBUG_NO_INLINE
#include "EmbeddedBoundaryShapeSphere.I"
#endif

namespace SAMRAI {
   namespace appu {

template<int DIM>
EmbeddedBoundaryShapeSphere<DIM>::EmbeddedBoundaryShapeSphere(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = object_name;

   tbox::IEEE::initializeArrayToSignalingNaN(d_center,DIM);
   tbox::IEEE::setNaN(d_radius);

   getFromInput(input_db);
}

template<int DIM>
EmbeddedBoundaryShapeSphere<DIM>::~EmbeddedBoundaryShapeSphere()
{  
}

template<int DIM> void
EmbeddedBoundaryShapeSphere<DIM>::printClassData(ostream& os) const
{
   os << "d_object_name = " << d_object_name << endl;
   os << "d_radius = " << d_radius << endl;
   for (int i = 0; i < DIM; i++) {
      os << "d_center[" << i << "] = " << d_center[i] << endl;
   }
   
}

template<int DIM> 
void EmbeddedBoundaryShapeSphere<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   /*
    * MUST supply a center and radius.
    */
   d_radius = db->getDouble("radius");
   
   tbox::Array<double> temp_center;
   temp_center = db->getDoubleArray("center");
   for (int i=0; i < DIM; i++) {
      d_center[i] = temp_center[i];
   }
}



}
}
#endif
