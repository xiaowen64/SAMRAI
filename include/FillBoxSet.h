//
// File:	FillBoxSet.h
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Utility class for "smart" boxlist operations in comm schedules.
//

#ifndef included_xfer_FillBoxSet
#define included_xfer_FillBoxSet

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_iostream
#include <iostream>
using namespace std;
#define included_iostream
#endif

namespace SAMRAI {
   namespace xfer {

/**
 * Class FillBoxSet is a utility class that provides a "smart" 
 * set of boxlist operations for communication schedules.  Specifically, 
 * it provides some of the basic box calculus operations, such as set 
 * intersection and difference, found in the simple boxlist class.  However, 
 * this class is used in very specific situations involving "fill boxes"
 * in communication schedule creation, to avoid potentially expensive 
 * exhaustive searches.  Essentially, it limits searches for box 
 * intersections to those that intersect the bounding box.
 *
 * @see hier::Box
 * @see hier::BoxList
 */

template<int DIM> 
class FillBoxSet
{
public:
   /**
    * The constructor defines a bounding box for the fill box set to 
    * be the argument box and creates a box set with that one box in it.  
    */
   FillBoxSet(const hier::Box<DIM>& box);

   /**
    * Create a new fill box set and copy the bounding box and
    * boxes from the argument fill box set.
    */
   FillBoxSet(const FillBoxSet<DIM>& fill_box_set);

   /**
    * Default constructor creates a new fill box set with an empty 
    * bounding box and an empty box set.  The box information must be 
    * set by calling one of the resetFillBoxes() functions to use the 
    * constructed object.
    */
   FillBoxSet();

   /**
    * The destructor releases all fill box set storage.
    */
   virtual ~FillBoxSet();

   /**
    * Copy bounding box and boxes from the argument fill box set.
    */
   virtual FillBoxSet& operator=(const FillBoxSet& fill_box_set);

   /**
    * Reset bounding box to be the argument box and replace fill box set 
    * with a set containing that one box.
    */
   virtual void resetFillBoxes(const hier::Box<DIM>& box);

   /**
    * Reset fill box set to be the argument box list and set bounding box 
    * to be the bounding box of the list.
    */
   virtual void resetFillBoxes(const hier::BoxList<DIM>& boxes);

   /**
    * Add new fill box to head of box list representing the fill box set.
    */
   void addFillBox(const hier::Box<DIM>& box);

   /**
    * Return const reference to bounding box for all boxes in the fill box set.
    */
   const hier::Box<DIM>& getBoundingBox();

   /**
    * Return const reference to box list describing the fill box set.
    */
   const hier::BoxList<DIM>& getBoxList() const;

   /**
    * Return (non-const) reference to box list describing the fill box set. 
    * This function should be used with care as the box list can be changed. 
    */
   hier::BoxList<DIM>& getBoxListToChange();

   /**
    * Return number of boxes in fill box set.
    */
   int getNumberOfBoxes() const;

   /**
    * Invoke box list routine of same name on the fill box set owned by
    * this object.
    */
   void simplifyBoxes();

   /**
    * If argument box intersects the bounding box owned by this object, 
    * then invoke the box list routine of same name on the fill box set 
    * owned by this object.  Otherwise do nothing.
    */
   void removeIntersections(const hier::Box<DIM>& takeaway);

   /**
    * Invoke box list routine of same name on the fill box set owned by
    * this object.
    */
   void removeIntersections(const hier::BoxList<DIM>& takeaway);

   /**
    * If argument box intersects the bounding box, then invoke the box 
    * list routine of same name on the fill box set owned by this object.  
    * Otherwise do nothing.
    */
   void intersectBoxes(const hier::Box<DIM>& box);

   /**
    * Invoke box list routine of same name on the fill box set owned by
    * this object.
    */
   void intersectBoxes(const hier::BoxList<DIM>& boxes);

   /**
    * Print all class member data for this fill box set object
    * to specified output stream.
    */
   virtual void print(ostream& os = tbox::plog) const;

private:
   hier::Box<DIM> d_bounding_box;
   hier::BoxList<DIM> d_boxes;
   bool d_recompute_bounding_box;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "FillBoxSet.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FillBoxSet.C"
#endif
