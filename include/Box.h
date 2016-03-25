//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/Box.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Box representing a portion of the AMR index space
//

#ifndef included_hier_Box
#define included_hier_Box

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "Index.h"
#include "IntVector.h"
#include "tbox/DatabaseBox.h"

namespace SAMRAI {
   namespace hier {

template<int DIM> class BoxIterator;

/**
 * Class Box represents a n-dimensional box in the AMR index
 * space.  It is defined by lower and upper bounds given by index objects.
 * The box semantics assumes that the box is cell-centered.  A cell-centered
 * convention implies that the index set covered by the box includes both
 * the lower and upper bounds.
 *
 * Class Box<DIM> is translated into classes Box1, Box2, and
 * Box3 after being passed through a preprocessor.
 *
 * @see hier::BoxIterator 
 * @see hier::Index
 */

template<int DIM> class Box
{
public:
   /**
    * The default constructor creates an ``empty'' box.
    */
   Box();

   /**
    * Create a box describing the index space between lower and upper.  The
    * box is assumed to be cell centered and include all elements between lower
    * and upper, including the end points.
    */
   Box(const Index<DIM>& lower, const Index<DIM>& upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   Box(const Box<DIM>& box);

   /**
    * Construct a Box<DIM> from a DatabaseBox.
    */
   Box(const tbox::DatabaseBox& box);
   
   /**
    * The destructor for Box.
    */
   ~Box();

   /**
    * The assignment operator copies the index space of the argument box.
    */
   Box<DIM>& operator=(const Box<DIM>& box);

   /**
    * Return a non-const lower index of the box.
    */
   Index<DIM>& lower();

   /**
    * Return a non-const upper index of the box.
    */
   Index<DIM>& upper();

   /**
    * Return a const lower index of the box.
    */
   const Index<DIM>& lower() const;

   /**
    * Return a const upper index of the box.
    */
   const Index<DIM>& upper() const;

   /** 
    * Return the i'th component (non-const) of the lower index of the box.
    */
   int& lower(const int i);

   /** 
    * Return the i'th component (non-const) of the upper index of the box.
    */
   int& upper(const int i);

   /** 
    * Return the i'th component (const) of the lower index of the box.
    */
   const int& lower(const int i) const;

   /** 
    * Return the i'th component (const) of the upper index of the box.
    */
   const int& upper(const int i) const;

   /**
    * Return whether the box is ``empty''.  A box is empty if any of the
    * lower bounds is greater than the corresponding upper bound.  An empty
    * box has a size of zero.
    */
   bool empty() const;

   /**
    * Set the index space represented by the box to empty.
    */
   void setEmpty();

   /**
    * Return the number of cells (an integer) represented by the box in
    * the given coordinate direction.
    */
   int numberCells(const int i) const;
 
   /**
    * Return the number of cells (a vector of integers) represented by
    * the box in every coordinate direction.
    */
   IntVector<DIM> numberCells() const;

   /**
    * Calculate the number of indices represented by the box.  If the box
    * is empty, then the number of index points within the box is zero.
    */
   int size() const;

   /**
    *  Return the dimension of the box that is longest.
    */
   int longestDimension() const;

   /**
    * Given an index, calculate the offset of the index into the box.
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for later array indexing operations.
    */
   int offset(const Index<DIM>& p) const;

   /**
    * Check whether an index lies within the bounds of the box.
    */
   bool contains(const Index<DIM>& p) const;

   /**
    * Check whether a given box lies within the bounds of the box.
    */
   bool contains(const Box<DIM>& b) const;

   /**
    * Check whether two boxes represent the same portion of index space.
    */
   int operator==(const Box<DIM>& box) const;

   /**
    * Check whether two boxes cover different portions of index space.
    */
   int operator!=(const Box<DIM>& box) const;

   /**
    * Calculate the intersection of the index spaces of two boxes.  The
    * intersection with an empty box always yields an empty box.
    */
   Box<DIM> operator*(const Box<DIM>& box) const;

   /**
    * Return true if two boxes have a non-empty intersection.  
    * Otherwise, return false. 
    */
   bool intersects(const Box<DIM>& box) const;

   /**
    * Calculate the bounding box of two boxes.  Note that this is not
    * the union of the two boxes (since union is not closed over boxes),
    * but rather the smallest box that contains both boxes.
    */
   Box<DIM> operator+(const Box<DIM>& box) const;

   /**
    * Increase the bounding box to include the argument box.
    */
   Box<DIM>& operator+=(const Box<DIM>& box);

   /**
    * Return true if this box can be coalesced with the argument box,
    * and set this box to the union of the boxes.  Otherwise, return false
    * and leave boxes as is.  Two boxes may be coalesced if their
    * union is a box (recall that index set union is not closed over boxes).  
    * If either box is empty, then the return value is true and this box 
    * becomes the union of the two.
    */
   bool coalesceWith(const Box<DIM>& box);

   /**
    * Grow a box by the specified ghost cell width.  The lower bound is
    * decremented by the width, and the upper bound is incremented by the
    * width.  All dimensions are grown by the corresponding component in
    * the IntVector; ghost cell widths may be different in each dimension.
    * Negative ghost cell widths will shrink the box.
    */
   void grow(const IntVector<DIM>& ghosts);

   /**
    * Grow a box by the specified ghost cell width in the given coordinate
    * direction in index space.  The lower bound is decremented by the 
    * width, and the upper bound is incremented by the width.  Note that
    * negative ghost cell widths will shrink the box.
    */
   void grow(const int direction, const int ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * directions (i.e., only lower index is changed).
    */
   void growLower(const IntVector<DIM>& ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * bound of given direction in index space.
    */
   void growLower(const int direction, const int ghosts);
  
   /**
    * Similar to grow() function. However, box is only grown in upper
    * directions (i.e., only upper index is changed).
    */
   void growUpper(const IntVector<DIM>& ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in upper
    * bound of given direction in index space.
    */
   void growUpper(const int direction, const int ghosts);

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * lengthened (never shortened).  The sign of @c ghosts refer to whether
    * the box is lengthened in the upper or lower side.
    */
   void lengthen(const int direction, const int ghosts);

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * shortened (never lengthened).  The sign of @c ghosts refer to whether
    * the box is shortened in the upper or lower side.
    */
   void shorten(const int direction, const int ghosts);

   /**
    * Shift a box by the specified amount (a vector of integers).
    * The new box is located at (lower+offset, upper+offset).
    */
   void shift(const IntVector<DIM>& offset);

   /**
    * Similar to shift() function above, but shift occurs only in specified
    * direction in index space.  The new box is located at (lower+offset, 
    * upper+offset) in that direction.
    */
   void shift(const int direction, const int offset);

   /**
    * Rotate 90 degrees around origin.  Currently works for 2D only.
    */ 
   void rotate(const int rotation_number);

   /**
    * Refine the index space of a box by specified vector ratio.  Each
    * component of the box is multiplied by the refinement ratio,
    * then @c (ratio-1) is added to the upper corner.
    */
   void refine(const IntVector<DIM>& ratio);

   /**
    * Coarsen the index space of a box by specified vector ratio.  Each
    * component is divided by the specified coarsening ratio and rounded
    * (if necessary) such that the coarsened box contains the cells that
    * are the parents of the refined box.  In other words, refining a
    * coarsened box will always yield a box that is equal to or larger
    * than the original box.
    */
   void coarsen(const IntVector<DIM>& ratio);

   /**
    * This assignment operator constructs a Box<DIM> given a DatabaseBox.
    */
   Box<DIM>& operator=(const tbox::DatabaseBox& box);

   /**
    * Sets a Box<DIM> from a tbox::DatabaseBox and returns a reference to
    * the Box<DIM>.
    */
   Box<DIM>& Box_from_DatabaseBox(const tbox::DatabaseBox& box);

   /**
    * Sets a Box<DIM> from a DatabaseBox. 
    */
   void set_Box_from_DatabaseBox(const tbox::DatabaseBox& box);

   /**
    * Returns a tbox::DatabaseBox generated from a Box<DIM>.
    */
   tbox::DatabaseBox DatabaseBox_from_Box() const;
   
   /**
    * Type conversion from Box<DIM> to Box
    */
   operator tbox::DatabaseBox();

   /**
    * Type conversion from Box<DIM> to Box
    */
   operator tbox::DatabaseBox() const;

   /**
    * Utility function to grow a box by the specified vector ghost cell
    * width.  A new box is returned and the argument is not changed.
    */
   static Box<DIM> grow(const Box<DIM>& box, const IntVector<DIM>& ghosts);

   /**
    * Utility function to shift a box by the specified offset.  A new
    * box is returned and the argument is not changed.
    */
   static Box<DIM> shift(const Box<DIM>& box, const IntVector<DIM>& offset);

   /**
    * Utility function to refine the index space of a box by the specified
    * refinement ratio.  A new box is returned and the argument is not changed.
    */
   static Box<DIM> refine(const Box<DIM>& box, const IntVector<DIM>& ratio);

   /**
    * Utility function to coarsen the index space of a box by the specified
    * coarsening ratio.  A new box is returned and the argument is not changed.
    */
   static Box<DIM> coarsen(const Box<DIM>& box, const IntVector<DIM>& ratio);
 
   /**
    * Read the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   template <int D> 
      friend std::istream& operator >> (std::istream& s, Box<D>& box);

   /**
    * Output the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   template <int D> 
      friend std::ostream& operator << (std::ostream& s, const Box<D>& box);

   /**
    * A box iterator iterates over the elements of a box.  This class is
    * defined elsewhere, and the typedef is used to point to that class.
    */
   typedef BoxIterator<DIM> Iterator;

private:
   static int coarsen(const int index, const int ratio);

   static bool coalesceIntervals(const int* lo1, const int* hi1,
			  const int* lo2, const int* hi2,
			  const int dim);

   void rotateAboutAxis(const int axis, const int num_rotations);

   Index<DIM> d_lo;
   Index<DIM> d_hi;
};

/**
 * Class BoxIterator is an iterator that provides methods for
 * stepping through the index space associated with a box.  The indices
 * are enumerated in column-major (e.g., Fortran) order.  The iterator
 * should be used as follows:
   \verbatim
   Box<DIM> box;
   ...
   for (Box<DIM>::Iterator b(box); b; b++) {
      // use index b of the box
   }
   \endverbatim
 * Note that the box iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see hier::Index
 * @see hier::Box
 */

template<int DIM> class BoxIterator

{
public:
   /**
    * Default constructor for the box iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    *
    * @see initialize().
    */
   BoxIterator();

   /**
    * Constructor for the box iterator.  The iterator will enumerate the
    * indices in the argument box.
    */
   BoxIterator(const Box<DIM>& box);

   /**
    * Copy constructor for the box iterator.
    */
   BoxIterator(const BoxIterator<DIM>& iterator);

   /**
    * Assignment operator for the box iterator.
    */
   BoxIterator<DIM>& operator=(const BoxIterator<DIM>& iterator);

   /**
    * Destructor for the box iterator.
    */
   ~BoxIterator<DIM>();

   /**
    * @brief Initializer for the box iterator.
    *
    * The iterator will enumerate the indices in the box.
    */
   void initialize(const Box<DIM>& box);

   /**
    * Return the current index in the box.  This operation is undefined
    * if the iterator is past the last Index in the box.
    */
   const Index<DIM>& operator*() const;

   /**
    * Return the current index in the box.  This operation is undefined
    * if the iterator is past the last Index in the box.
    */
   const Index<DIM>& operator()() const;

   /**
    * Return true if the iterator points to a valid index in the box.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION 
   /**
    * Return a non-NULL if the iterator points to a valid index in the box.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid item in the box.  This
    * operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next index in the box.
    */
   void operator++(int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool operator==(const BoxIterator<DIM>& iterator) const;
  
   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const BoxIterator<DIM>& iterator) const;

private:
   Index<DIM> d_index;
   Box<DIM> d_box;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "Box.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "Box.C"
#endif

#endif

