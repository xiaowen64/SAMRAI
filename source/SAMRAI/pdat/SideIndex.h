/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_SideIndex
#define included_pdat_SideIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class SideIndex implements a simple n-dimensional integer
 * vector for side centered variables.  Side indices contain an integer
 * index location in AMR index space along with the designated side axis
 * (X=0, Y=1, or Z=2).  See the side box geometry class for more information
 * about the mapping between the AMR index space and the side indices.
 *
 * @see hier::Index
 * @see pdat::SideData
 * @see pdat::SideGeometry
 * @see pdat::SideIterator
 */

class SideIndex:public hier::Index
{
public:
   /**
    * The default constructor for a side index creates an uninitialized index.
    */
   explicit SideIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a side index from a regular index, axis, and side.  The axis
    * can be one of SideIndex::X (0), SideIndex::Y (1), or
    * SideIndex::Z (2). The side argument can be one of the constants
    * SideIndex::Lower (0) or SideIndex::Upper (1).
    */
   explicit SideIndex(
      const hier::Index& rhs,
      const int axis,
      const int side);

   /**
    * The copy constructor creates a side index equal to the argument.
    */
   SideIndex(
      const SideIndex& rhs);

   /**
    * The assignment operator sets the side index equal to the argument.
    */
   SideIndex&
   operator = (
      const SideIndex& rhs);

   /**
    * The side index destructor does nothing interesting.
    */
   ~SideIndex();

   /**
    * Get the axis for which this side index is defined (X=0, Y=1, Z=2).
    */
   int
   getAxis() const;

   /**
    * Set the side axis (X=0, Y=1, Z=2).
    */
   void
   setAxis(
      const int axis);

   /**
    * Convert the side index into the index on the left hand side
    * (argument side == 0) or the right hand side (argument side == 1).
    */
   hier::Index
   toCell(
      const int side) const;

   /**
    * Plus-equals operator for a side index and an integer vector.
    */
   SideIndex&
   operator += (
      const hier::IntVector& rhs);

   /**
    * Plus operator for a side index and an integer vector.
    */
   SideIndex
   operator + (
      const hier::IntVector& rhs) const;

   /**
    * Plus-equals operator for a side index and an integer.
    */
   SideIndex&
   operator += (
      const int rhs);

   /**
    * Plus operator for a side index and an integer.
    */
   SideIndex
   operator + (
      const int rhs) const;

   /**
    * Minus-equals operator for a side index and an integer vector.
    */
   SideIndex&
   operator -= (
      const hier::IntVector& rhs);

   /**
    * Minus operator for a side index and an integer vector.
    */
   SideIndex
   operator - (
      const hier::IntVector& rhs) const;

   /**
    * Minus-equals operator for a side index and an integer.
    */
   SideIndex&
   operator -= (
      const int rhs);

   /**
    * Minus operator for a side index and an integer.
    */
   SideIndex
   operator - (
      const int rhs) const;

   /**
    * Times-equals operator for a side index and an integer vector.
    */
   SideIndex&
   operator *= (
      const hier::IntVector& rhs);

   /**
    * Times operator for a side index and an integer vector.
    */
   SideIndex
   operator * (
      const hier::IntVector& rhs) const;

   /**
    * Times-equals operator for a side index and an integer.
    */
   SideIndex&
   operator *= (
      const int rhs);

   /**
    * Times operator for a side index and an integer.
    */
   SideIndex
   operator * (
      const int rhs) const;

   /**
    * Returns true if two side index objects are equal.  All components
    * and the corresponding side axes must be the same for equality.
    */
   bool
   operator == (
      const SideIndex& rhs) const;

   /**
    * Returns true if two side index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    */
   bool
   operator != (
      const SideIndex& rhs) const;

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
#include "SAMRAI/pdat/SideIndex.I"
#endif
#endif
