//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/embedded_boundary/EmbeddedBoundaryShapeSphere.h $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Sphere embedded boundary shape
//              
// 

#ifndef included_appu_EmbeddedBoundaryShapeSphere
#define included_appu_EmbeddedBoundaryShapeSphere

#include "SAMRAI_config.h"

#include "tbox/Database.h"
#include "EmbeddedBoundaryShape.h"
#include "EmbeddedBoundaryDefines.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

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

 * \endverbatim
 *
 */
      
template<int DIM>
class EmbeddedBoundaryShapeSphere : public EmbeddedBoundaryShape<DIM>
{
public:
   
   /*!
    * The constructor initializes center and radius to NaN.
    *
    * @param object_name name of object of this class
    * @param input_db    the input database which contains radius and 
    *                    center specification.
    */
   EmbeddedBoundaryShapeSphere(const std::string& object_name,
                               tbox::Pointer<tbox::Database> input_db);
   
   /*!
    * The destructor does nothing.
    */
   ~EmbeddedBoundaryShapeSphere<DIM>();

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShape base class.  This method indicates 
    * whether the supplied xyz coordinates are inside or outside of 
    * the sphere.
    *
    * @param xyz  double array[DIM] specifying coordinates. 
    */
   bool isInside(const double* xyz) const;

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShape base class.  This method indicates 
    * whether the array of xyz coordinates are inside or outside of 
    * the sphere.
    *
    * @param nx integer array [DIM] specifying number of points in each dir
    * @param dx double array [DIM] specifying spacing of points in each dir
    * @param origin double array [DIM] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points 
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   void isInside(const int* nx,
                 const double* dx,
                 const double* origin,
                 int* inout) const;
   
   /*!
    * Dump data to supplied stream.
    */
   virtual void printClassData(std::ostream& os) const;
   
private:   

   /*
    * Read name, center, and radius information from input.  The name
    * is optional but center and radius MUST be specified in the input 
    * file.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   std::string d_object_name;

   /*
    * Center and radius of the sphere.
    */
   double d_center[DIM];
   double d_radius;

};   
 
     
}
}

#ifndef DEBUG_NO_INLINE
#include "EmbeddedBoundaryShapeSphere.I"
#endif

#endif // included_EmbeddedBoundaryShapeSphere

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EmbeddedBoundaryShapeSphere.C"
#endif
   
