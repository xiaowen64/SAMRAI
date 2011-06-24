/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Sphere embedded boundary shape 
 *
 ************************************************************************/

#ifndef included_appu_EmbeddedBoundaryShapeSphere_C
#define included_appu_EmbeddedBoundaryShapeSphere_C

#include "SAMRAI/appu/EmbeddedBoundaryShapeSphere.h"

#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/appu/EmbeddedBoundaryShapeSphere.I"
#endif

namespace SAMRAI {
namespace appu {

EmbeddedBoundaryShapeSphere::EmbeddedBoundaryShapeSphere(
   const tbox::Dimension& dim,
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db):
   d_dim(dim),
   d_object_name(object_name),
   d_radius(tbox::MathUtilities<double>::getSignalingNaN())
{
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, d_dim.getValue());

   getFromInput(input_db);
}

EmbeddedBoundaryShapeSphere::~EmbeddedBoundaryShapeSphere()
{
}

void
EmbeddedBoundaryShapeSphere::printClassData(
   std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_radius = " << d_radius << std::endl;
   for (int i = 0; i < d_dim.getValue(); i++) {
      os << "d_center[" << i << "] = " << d_center[i] << std::endl;
   }

}

void EmbeddedBoundaryShapeSphere::getFromInput(
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
   for (int i = 0; i < d_dim.getValue(); i++) {
      d_center[i] = temp_center[i];
   }
}

}
}
#endif
