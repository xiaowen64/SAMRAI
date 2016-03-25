//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/face/FaceIndex.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Index for face centered patch data types
//

#ifndef included_pdat_FaceIndex
#define included_pdat_FaceIndex

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "Index.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class FaceIndex<DIM> implements a simple n-dimensional integer
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

template<int DIM> class FaceIndex : public hier::Index<DIM>
{
public:
   /**
    * The default constructor for a face index creates an uninitialized index.
    */
   FaceIndex();

   /**
    * Construct a face index from a regular index, axis, and face.  The axis
    * can be one of FaceIndex<DIM>::X (0), FaceIndex<DIM>::Y (1), or
    * FaceIndex<DIM>::Z (2). The face argument can be one of the constants
    * FaceIndex<DIM>::Lower (0) or FaceIndex<DIM>::Upper (1).
    */
   FaceIndex(const hier::Index<DIM>& rhs, const int axis, const int face);

   /**
    * The copy constructor creates a face index equal to the argument.
    */
   FaceIndex(const FaceIndex<DIM>& rhs);

   /**
    * The assignment operator sets the face index equal to the argument.
    */
   FaceIndex<DIM>& operator=(const FaceIndex<DIM>& rhs);

   /**
    * The face index destructor does nothing interesting.
    */
   ~FaceIndex<DIM>();

   /**
    * Get the axis for which this face index is defined (X=0, Y=1, Z=2).
    */
   int getAxis() const;

   /**
    * Set the face axis (X=0, Y=1, Z=2).
    */
   void setAxis(const int axis);

   /**
    * Convert the face index into the index on the left hand face 
    * (argument face == 0) or the right hand face (argument face == 1).
    */
   hier::Index<DIM> toCell(const int face) const;

   /**
    * Plus-equals operator for a face index and an integer vector.
    */
   FaceIndex<DIM>& operator+=(const hier::IntVector<DIM>& rhs);

   /**
    * Plus operator for a face index and an integer vector.
    */
   FaceIndex<DIM> operator+(const hier::IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for a face index and an integer.
    */
   FaceIndex<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for a face index and an integer.
    */
   FaceIndex<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for a face index and an integer vector.
    */
   FaceIndex<DIM>& operator-=(const hier::IntVector<DIM>& rhs);

   /**
    * Minus operator for a face index and an integer vector.
    */
   FaceIndex<DIM> operator-(const hier::IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for a face index and an integer.
    */
   FaceIndex<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for a face index and an integer.
    */
   FaceIndex<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for a face index and an integer vector.
    */
   FaceIndex<DIM>& operator*=(const hier::IntVector<DIM>& rhs);

   /**
    * Times operator for a face index and an integer vector.
    */
   FaceIndex<DIM> operator*(const hier::IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for a face index and an integer.
    */
   FaceIndex<DIM>& operator*=(const int rhs);

   /**
    * Times operator for a face index and an integer.
    */
   FaceIndex<DIM> operator*(const int rhs) const;

   /**
    * Returns true if two face index objects are equal.  All components
    * and the corresponding face axes must be the same for equality.
    */
   bool operator==(const FaceIndex<DIM>& rhs) const;

   /**
    * Returns true if two face index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    */
   bool operator!=(const FaceIndex<DIM>& rhs) const;

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
#ifndef DEBUG_NO_INLINE
#include "FaceIndex.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceIndex.C"
#endif
