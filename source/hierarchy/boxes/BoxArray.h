//
// File:	BoxArray.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An array of boxes that complements BoxList
//

#ifndef included_hier_BoxArray
#define included_hier_BoxArray

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_tbox_DatabaseBox
#include "tbox/DatabaseBox.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> class BoxList;

/**
 * Class BoxArray is a container represents an array of boxes.
 * Once the array is created, it can be resized.  However, it is usually
 * better to use the BoxList class to represent collections of boxes 
 * that change in number.
 *
 * @see hier::Box
 * @see hier::BoxList
 */

template<int DIM> class BoxArray
{
public:
   /**
    * Create an array of boxes with space for n boxes.  All boxes are
    * initialized to empty.
    */
   BoxArray(const int n = 0);

   /**
    * Create a box array and copy box data from the array argument.
    */
   BoxArray(const tbox::Array< Box<DIM> >& array);

   /**
    * The const constructor creates an array of boxes and copies
    * box data from the argument box array.
    */
   BoxArray(const BoxArray<DIM>& array);

   /**
    * Create a box array and copy box data from the list.
    */
   BoxArray(const BoxList<DIM>& list);

   /**
    * Create a regular box array from an array of tbox::DatabaseBox objects.
    */
   BoxArray(const tbox::Array<tbox::DatabaseBox>& array);

   /**
    * Type conversion from BoxArray<DIM> to tbox::Array<tbox::DatabaseBox>.
    */
   operator tbox::Array<tbox::DatabaseBox>() const;
   
   /**
    * Create a box array using boxes in tbox::Array<tbox::DatabaseBox> for the data.
    */
   BoxArray<DIM>& operator=(const tbox::Array<tbox::DatabaseBox>& array);

   /**
    * Create a box array and copy data from the array argument.
    */
   BoxArray<DIM>& operator=(const BoxArray<DIM>& array);

   /**
    * Create a box array and copy data from the list argument.
    */
   BoxArray<DIM>& operator=(const BoxList<DIM>& list);

   /**
    * The BoxArray destructor releases the box array data.
    */
   ~BoxArray<DIM>();

   /**
    * Return the number of boxes in the array.
    */
   int getNumberOfBoxes() const;

   /**
    * Return the number of boxes in the array.  Identical to getNumberOfBoxes(),
    * but this method is common to several container classes.
    */
   int size() const;

   /**
    * Return a reference to the i-th box.  No bounds checking.
    */
   Box<DIM>& getBox(const int i);

   /**
    * Return a const reference to the i-th box.  No bounds checking.
    */
   const Box<DIM>& getBox(const int i) const;

   /**
    * Return a reference to the i-th box.  No bounds checking.
    */
   Box<DIM>& operator()(const int i);

   /**
    * Return a const reference to the i-th box.  No bounds checking.
    */
   const Box<DIM>& operator()(const int i) const;

   /**
    * Create a BoxArray<DIM> from a tbox::Array<tbox::DatabaseBox>.
    */
   BoxArray<DIM>& BoxArray_from_Array(tbox::Array<tbox::DatabaseBox> array);

   /**
    * Sets a BoxArray<DIM> from a tbox::Array<tbox::DatabaseBox>.
    */
   void set_BoxArray_from_Array(tbox::Array<tbox::DatabaseBox> array);

   /**
    * Check whether an index lies within the bounds of the collection
    * of boxes.
    */
   bool contains(const Index<DIM>& p) const;

   /**
    * Grow all boxes in the box array by the specified ghost cell width.
    */
   void grow(const IntVector<DIM>& ghosts);

   /**
    * Shift all boxes in the box array by the specified offset.
    */
   void shift(const IntVector<DIM>& offset);

   /**
    * rotate the box array clockwise 90 degrees times the rotation number.
    * Currently works only in 2D.
    */
   void rotate(int rotation_number);

   /**
    * Refine the index space of each box in the box array by the
    * specified vector refinement ratio.
    */
   void refine(const IntVector<DIM>& ratio);

   /**
    * Coarsen the index space of each box in the box array by the
    * specified vector coarsening ratio.
    */
   void coarsen(const IntVector<DIM>& ratio);

   /**
    * Resize the array to contain the specified number of boxes.
    * The resizing operation copies the old box values into the
    * new boxes where the indices in the array overlap.
    */
   void resizeBoxArray(const int n);

   /**
    * Print all boxes in this array to specified output stream.
    */
   void print(ostream& os = tbox::plog) const;

private:
   tbox::Array< Box<DIM> > d_boxes;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "BoxArray.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxArray.C"
#endif

#endif

