//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/embedded_boundary/EmbeddedBoundaryShapeSphere.C $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Sphere embedded boundary shape
//              
// 

#ifndef included_appu_EmbeddedBoundaryShapeSphere_C
#define included_appu_EmbeddedBoundaryShapeSphere_C

#include "EmbeddedBoundaryShapeSphere.h"


#include "tbox/MathUtilities.h"

#ifdef DEBUG_NO_INLINE
#include "EmbeddedBoundaryShapeSphere.I"
#endif

namespace SAMRAI {
   namespace appu {

template<int DIM>
EmbeddedBoundaryShapeSphere<DIM>::EmbeddedBoundaryShapeSphere(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = object_name;

   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, DIM);
   d_radius = tbox::MathUtilities<double>::getSignalingNaN();

   getFromInput(input_db);
}

template<int DIM>
EmbeddedBoundaryShapeSphere<DIM>::~EmbeddedBoundaryShapeSphere()
{  
}

template<int DIM> void
EmbeddedBoundaryShapeSphere<DIM>::printClassData(std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_radius = " << d_radius << std::endl;
   for (int i = 0; i < DIM; i++) {
      os << "d_center[" << i << "] = " << d_center[i] << std::endl;
   }
   
}

template<int DIM> 
void EmbeddedBoundaryShapeSphere<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
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
