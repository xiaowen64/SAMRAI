/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for analytic embedded boundaries 
 *
 ************************************************************************/

#ifndef included_appu_EmbeddedBoundaryShape
#define included_appu_EmbeddedBoundaryShape

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Array.h"

namespace SAMRAI {
namespace appu {

/*!
 * @brief An abstract base class from which the different embedded boundary
 * analytic shapes used in SAMRAI are derived.  It specifies virtual
 * implementations of several functions used to define an analytic
 * embedded boundary shape.
 */
class EmbeddedBoundaryShape
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
    * @param xyz double array dimension [tbox::Dimension::MAXIMUM_DIMENSION_VALUE], specifying the location of
    *            a single point.
    */
   virtual bool
   isInside(
      const double* xyz) const;

   /*!
    * Same implementation as above, but instead takes an array of points.
    * Unlike the method above, which returns a boolean for a single point,
    * this method sets an integer for any array of points layed out in a
    * patch.  The "nx" argument specifies the dimension of the array of
    * points in DIM directions - it should be an int dimensioned
    * [tbox::Dimension::MAXIMUM_DIMENSION_VALUE].  The "dx" argument specifies the grid spacing in each
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
    * @param nx integer array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying number of points in each dir
    * @param dx double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying spacing of points in each dir
    * @param origin double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   virtual void
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
      std::ostream& os) const = 0;

private:
};

}
}
#endif
