//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/embedded_boundary/EmbeddedBoundaryShapePolygon.h $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Polygon embedded boundary shape
//              
// 

#ifndef included_appu_EmbeddedBoundaryShapePolygon
#define included_appu_EmbeddedBoundaryShapePolygon

#include "SAMRAI_config.h"

#include "tbox/Database.h"
#include "EmbeddedBoundaryDefines.h"
#include "EmbeddedBoundaryShape.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
   namespace appu {
   
/*!
 * @brief Provides an analytic description of a convex polygon.
 * It inherets from the EmbeddedBoundaryShape base class and provides a 
 * concrete implementation of the "isInside()" method, which specifies 
 * whether a cell is INSIDE the convex poly.
 *
 * The user must specify at least three coordinates that define the vertices
 * the poly.  If the problem is 3D, a height must also be specified.  An
 * example input entry would look like:
 *
 * \verbatim
 *        Polygon1{
 *           type = "POLYGON"
 *           coords_1 = 1.0 , 1.0
 *           coords_2 = 2.0 , 1.0
 *           coords_3 = 2.0 , 2.0
 *           coords_4 = 1.0 , 2.0
 *           height = 8.0
 *        }
 * \endverbatim
 *
 */

template<int DIM>
class EmbeddedBoundaryShapePolygon : public EmbeddedBoundaryShape<DIM>
{
public:
   /*!
    * @param object_name name of object of this class
    * @param input_db    the input database which contains radius and 
    *                    center specification.
    */
   EmbeddedBoundaryShapePolygon(const std::string& object_name,
                tbox::Pointer<tbox::Database> input_db);
   
   /*!
    * The destructor does nothing.
    */
   ~EmbeddedBoundaryShapePolygon();

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShapeX base class.  This method indicates 
    * whether the supplied xyz coordinates are inside or outside of 
    * the polygon.
    *
    * @param xyz  double array[DIM] specifying coordinates. 
    */
   bool isInside(const double* xyz) const;

   /*!
    * Concrete implementation of the isInside() method defined by the
    * EmbeddedBoundaryShapeX base class.  This method indicates 
    * whether the array of xyz coordinates are inside or outside of 
    * the polygon.
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

   /*!
    * Read name, and vertex information from input.  The name is optional but
    * at least three vertices must be specified by the input file.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   /*!
    * Returns TRUE if points p1 and p2 are on the same side of line segment ab;
    * FALSE otherwise.
    */
   bool sameSide(double p1[3], double p2[3],
                 double a[3], double b[3]) const;

   /*!
    * Returns TRUE if a point p is within a convex polygon defined by vertices
    * (v_x[i], * v_y[i]), i = 0..d_num_vertices-1; FALSE otherwise.
    */
   bool pointInPolygon(tbox::Array<double> v_x,
                       tbox::Array<double> v_y,
                       double p_x,
                       double p_y) const;

   /*!
    * This method computes the cross prodcut, A = B cross C. 
    */
   void crossProduct(double a[3],
                     const double b[3],
                     const double c[3]) const;
   

   /*!
    * This method computes the dot product of A and B. 
    */
   double dotProduct(const double a[3],
                     const double b[3]) const;
   

   std::string d_object_name;

   /*
    * Coordinates in X,Y space for the polygon, and height in Z.
    */
   double d_height;

   /*
    * Machine roundoff
    */
   double d_eps;

   /*
    * Arrays of x and y vertices respectively.
    */
   tbox::Array<double> d_vx;
   tbox::Array<double> d_vy;

   int d_num_vertices;
   
};   

} // namespace appu
} // namespace SAMRAI

#ifndef DEBUG_NO_INLINE
#include "EmbeddedBoundaryShapePolygon.I"
#endif

#endif // included_EmbeddedBoundaryShapePolygon

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EmbeddedBoundaryShapePolygon.C"
#endif

   
