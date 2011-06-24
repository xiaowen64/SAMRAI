/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Sphere embedded boundary shape 
 *
 ************************************************************************/

#ifndef included_appu_EmbeddedBoundaryShapeSphere
#define included_appu_EmbeddedBoundaryShapeSphere

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/EmbeddedBoundaryShape.h"
#include "SAMRAI/appu/EmbeddedBoundaryDefines.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace appu {

/*!
 * @brief Provides an analytic description of a sphere.  It inherets
 * from the EmbeddedBoundaryShape base class and provides a concrete
 * implementation of the "isInside()" method, which specifies whether a
 * cell is INSIDE the sphere.
 *
 * The user must specify in the input a "center" and a "radius". An
 * example input entry would look like:
 *
 * \verbatim
 *        Shape1 {
 *           type = "SPHERE"
 *           center = 40.0 , 15.0, 15.0
 *           radius = 5.0
 *        }
 *
 * \endverbatim
 *
 */

class EmbeddedBoundaryShapeSphere:public EmbeddedBoundaryShape
{
public:
   /*!
    * The constructor initializes center and radius to NaN.
    *
    * @param object_name name of object of this class
    * @param input_db    the input database which contains radius and
    *                    center specification.
    */
   EmbeddedBoundaryShapeSphere(
      const tbox::Dimension& dim,
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db);

   /*!
    * The destructor does nothing.
    */
   ~EmbeddedBoundaryShapeSphere();

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShape base class.  This method indicates
    * whether the supplied xyz coordinates are inside or outside of
    * the sphere.
    *
    * @param xyz  double array[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying coordinates.
    */
   bool
   isInside(
      const double* xyz) const;

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShape base class.  This method indicates
    * whether the array of xyz coordinates are inside or outside of
    * the sphere.
    *
    * @param nx integer array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying number of points in each dir
    * @param dx double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying spacing of points in each dir
    * @param origin double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   void
   isInside(
      const int* nx,
      const double* dx,
      const double* origin,
      int* inout) const;

   /*!
    * Dump data to supplied stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const;

private:
   /*
    * Read name, center, and radius information from input.  The name
    * is optional but center and radius MUST be specified in the input
    * file.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db);

   const tbox::Dimension& d_dim;

   std::string d_object_name;

   /*
    * Center and radius of the sphere.
    */
   double d_center[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double d_radius;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/appu/EmbeddedBoundaryShapeSphere.I"
#endif

#endif // included_EmbeddedBoundaryShapeSphere
