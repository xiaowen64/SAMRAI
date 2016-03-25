//
// File:        EmbeddedBoundaryShapeSphere.h
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 726 $
// Modified:    $Date: 2005-11-10 17:19:40 -0800 (Thu, 10 Nov 2005) $
// Description: Sphere embedded boundary shape
//              
// 

#ifndef included_appu_EmbeddedBoundaryShapeSphere
#define included_appu_EmbeddedBoundaryShapeSphere

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_appu_EmbeddedBoundaryShape
#include "EmbeddedBoundaryShape.h"
#endif
#ifndef included_appu_EmbeddedBoundaryDefines
#include "EmbeddedBoundaryDefines.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

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
   EmbeddedBoundaryShapeSphere(const string& object_name,
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
   void printClassData(ostream& os) const;
   
private:   

   /*
    * Read name, center, and radius information from input.  The name
    * is optional but center and radius MUST be specified in the input 
    * file.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   string d_object_name;

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
   
