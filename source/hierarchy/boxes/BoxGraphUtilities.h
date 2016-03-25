//
// File:        BoxGraphUtilities.h
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 601 $
// Modified:    $Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Utility class for operations that reduce complexity of box calculus
//

#ifndef included_hier_BoxGraphUtilities
#define included_hier_BoxGraphUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif


namespace SAMRAI {
   namespace hier {

/*!
 * Class BoxGraphUtilities is a utility class that provides
 * methods for "expanding" an array of boxes when some of the
 * boxes touch periodic boundaries.  If a box touches a periodic
 * boundary, then there are one or more "virtual" boxes that
 * "wrap around" to the opposite side(s) of the domain:
 * 
 * 
 * @verbatim
 * 
 *         ---------------------------
 *         |                         |
 *         |                         |
 *         |                         |
 *         |-------------            |-------------
 *         | this       |            | the virtual|
 *         | box touches|            | wrap around|
 *         | a periodic |            | box        |
 *         | boundary   |            |            |
 *         |-------------            |-------------
 *         |                         |
 *         |       domain            |
 *         |        box              |
 *         ---------------------------
 * 
 * @endverbatim
 * 
 * 
 * The BoxTop, BoxGraph, and possibly other classes require
 * as input arrays of boxes in which the "virtual" boxes are
 * explicitly represented.
 */

template<int DIM> struct BoxGraphUtilities
{
   /*!
    * @brief Returns an array of boxes that includes entries for each
    * box that touches a periodic boundary.
    * 
    * The shift array must either have zero length, or have the same length
    * as \b in_boxes, otherwise an unrecoverable error will be thrown.
    * 
    * @param out_boxes contains all items in the \b in_boxes array, and contains
    *                  one or more additional items for any box that touches
    *                  a periodic boundary.
    * @param in_boxes  array of input boxes.
    * @param shifts    shift information for each of the input boxes.
    */
   static void makeBoxesPlusPeriodicBoxes(
      BoxArray<DIM>& out_boxes,
      const BoxArray<DIM>& in_boxes,
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Returns an array of boxes that includes entries for each
    * box that touches a periodic boundary.
    * 
    * The shift array must either have zero length, or have the same length
    * as \b in_boxes, otherwise an unrecoverable error will be thrown.
    * 
    * 
    * @param out_boxes contains all items in the \b in_boxes array, and contains
    *                  one or more additional items for any box that touches
    *                  a periodic boundary.
    * @param out_indices contains an entry for each box in \b out_boxes;
    *                    the entry indicates the box, w.r.t \b in_boxes,
    *                    from which the box was derived.
    * @param in_boxes  array of input boxes.
    * @param shifts    shift information for each of the input boxes.
    */
   static void makeBoxesPlusPeriodicBoxes(
      BoxArray<DIM>& out_boxes,
      tbox::Array<int>& out_indices,
      const BoxArray<DIM>& in_boxes,
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Returns the sum of shifts[j].getNumberItems().
    * 
    * This function is called by makeBoxesPlusPeriodicBoxes().
    * 
    * @param shifts periodic shift information for each box.
    */
   static int countPeriodicBoxes(
      const tbox::Array< tbox::List< IntVector<DIM> > >& shifts);

   /*!
    * @brief Compare function for use with qsort when sorting integers
    * in ascending order.
    * 
    * Sample usage:
    * 
    * @verbatim
    *   intarray[len];
    *   ...
    *   qsort(array, len, sizeof(int), qsortIntCompare);
    * @endverbatim
    * 
    */
   static int qsortIntCompare(const void *v, const void *w);
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGraphUtilities.C"
#endif

#endif

