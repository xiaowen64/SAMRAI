//
// File:        EmbeddedBoundaryShape.h
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 726 $
// Modified:    $Date: 2005-11-10 17:19:40 -0800 (Thu, 10 Nov 2005) $
// Description: Base class for analytic embedded boundaries
//              
// 

#ifndef included_appu_EmbeddedBoundaryShape
#define included_appu_EmbeddedBoundaryShape

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif

namespace SAMRAI {
   namespace appu {

/*!
 * @brief An abstract base class from which the different embedded boundary 
 * analytic shapes used in SAMRAI are derived.  It specifies virtual
 * implementations of several functions used to define an analytic
 * embedded boundary shape.
 */
template<int DIM> class EmbeddedBoundaryShape
{
public:
   /*!
    * The constructor and destructor essentially do nothing.
    */
   EmbeddedBoundaryShape();
   
   virtual ~EmbeddedBoundaryShape();

   /*!
    * Virtual implementation of the isInside() method.  The class
    * that inherets from this base class must provide an implementation of
    * this function, specifying whether the supplied xyz coordinates are 
    * inside or outside of the particular shape being defined.
    *
    * @param xyz double array dimension [DIM], specifying the location of 
    *            a single point.
    */
   virtual bool isInside(const double* xyz) const;

   /*!
    * Same implementation as above, but instead takes an array of points.
    * Unlike the method above, which returns a boolean for a single point, 
    * this method sets an integer for any array of points layed out in a
    * patch.  The "nx" argument specifies the dimension of the array of
    * points in DIM directions - it should be an int dimensioned 
    * [DIM].  The "dx" argument specifies the grid spacing in each 
    * direction - it should be a double dimensioned DIM.  The "origin" 
    * argument specifies the origin, or lower corner of the patch - 
    * it also should be a double dimensioned DIM.  From these three
    * pieces of information, the coordinates of all points on the 
    * patch can be computed.
    *
    * The "inout" argument is the returned integer array which
    * will contain a definition of whether each point on the patch is 
    * inside (1) or outside (0) the geometry.
    *
    * @param nx integer array [DIM] specifying number of points in each dir
    * @param dx double array [DIM] specifying spacing of points in each dir
    * @param origin double array [DIM] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points 
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   virtual void isInside(const int* nx,
                         const double* dx,
                         const double* origin,
                         int* inout) const;

   /*!
    * Dump data to supplied stream.
    */
   virtual void printClassData(ostream& os) const = 0;

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EmbeddedBoundaryShape.C"
#endif
