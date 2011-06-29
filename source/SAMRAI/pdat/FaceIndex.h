/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_FaceIndex
#define included_pdat_FaceIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class FaceIndex implements a simple n-dimensional integer
 * vector for face centered variables.  Face indices contain an integer
 * index location in AMR index space along with the designated face axis
 * (X=0, Y=1, or Z=2).  See the face box geometry class for more information
 * about the mapping between the AMR index space and the face indices.
 *
 * @see hier::Index
 * @see pdat::FaceData
 * @see pdat::FaceGeometry
 * @see pdat::FaceIterator
 */

class FaceIndex:public hier::Index
{
public:
   /**
    * The default constructor for a face index creates an uninitialized index.
    */
   explicit FaceIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a face index from a regular index, axis, and face.  The axis
    * can be one of FaceIndex::X (0), FaceIndex::Y (1), or
    * FaceIndex::Z (2). The face argument can be one of the constants
    * FaceIndex::Lower (0) or FaceIndex::Upper (1).
    */
   explicit FaceIndex(
      const hier::Index& rhs,
      const int axis,
      const int face);

   /**
    * The copy constructor creates a face index equal to the argument.
    */
   FaceIndex(
      const FaceIndex& rhs);

   /**
    * The assignment operator sets the face index equal to the argument.
    */
   FaceIndex&
   operator = (
      const FaceIndex& rhs);

   /**
    * The face index destructor does nothing interesting.
    */
   ~FaceIndex();

   /**
    * Get the axis for which this face index is defined (X=0, Y=1, Z=2).
    */
   int
   getAxis() const;

   /**
    * Set the face axis (X=0, Y=1, Z=2).
    */
   void
   setAxis(
      const int axis);

   /**
    * Convert the face index into the index on the left hand face
    * (argument face == 0) or the right hand face (argument face == 1).
    */
   hier::Index
   toCell(
      const int face) const;

   /**
    * Plus-equals operator for a face index and an integer vector.
    */
   FaceIndex&
   operator += (
      const hier::IntVector& rhs);

   /**
    * Plus operator for a face index and an integer vector.
    */
   FaceIndex
   operator + (
      const hier::IntVector& rhs) const;

   /**
    * Plus-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator += (
      const int rhs);

   /**
    * Plus operator for a face index and an integer.
    */
   FaceIndex
   operator + (
      const int rhs) const;

   /**
    * Minus-equals operator for a face index and an integer vector.
    */
   FaceIndex&
   operator -= (
      const hier::IntVector& rhs);

   /**
    * Minus operator for a face index and an integer vector.
    */
   FaceIndex
   operator - (
      const hier::IntVector& rhs) const;

   /**
    * Minus-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator -= (
      const int rhs);

   /**
    * Minus operator for a face index and an integer.
    */
   FaceIndex
   operator - (
      const int rhs) const;

   /**
    * Times-equals operator for a face index and an integer vector.
    */
   FaceIndex&
   operator *= (
      const hier::IntVector& rhs);

   /**
    * Times operator for a face index and an integer vector.
    */
   FaceIndex
   operator * (
      const hier::IntVector& rhs) const;

   /**
    * Times-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator *= (
      const int rhs);

   /**
    * Times operator for a face index and an integer.
    */
   FaceIndex
   operator * (
      const int rhs) const;

   /**
    * Returns true if two face index objects are equal.  All components
    * and the corresponding face axes must be the same for equality.
    */
   bool
   operator == (
      const FaceIndex& rhs) const;

   /**
    * Returns true if two face index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    */
   bool
   operator != (
      const FaceIndex& rhs) const;

   enum {
      X = 0,
      Y = 1,
      Z = 2,
      Lower = 0,
      Upper = 1
   };

private:
   int d_axis;
};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/FaceIndex.I"
#endif
#endif
