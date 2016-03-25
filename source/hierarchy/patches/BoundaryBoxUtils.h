//
// File:	$URL$
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2006 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2350 $
// Modified:	$LastChangedDate: 2008-09-10 11:21:21 -0700 (Wed, 10 Sep 2008) $
// Description:	Generic utilities for boundary box calculus.
//

#ifndef included_hier_BoundaryBoxUtils
#define included_hier_BoundaryBoxUtils

#include "SAMRAI_config.h"
#include "BoundaryBox.h"


namespace SAMRAI {
   namespace hier {

/*!
 * @brief Perform shifts, extensions, etc on a BoundaryBox using the box's
 * location index and type.
 *
 * @see hier::BoundaryBox
 */

template<int DIM>
class BoundaryBoxUtils
{
public:

BoundaryBoxUtils();


/*!
 * @brief Construct with a boundary box.
 *
 * @see setBoundaryBox()
 */
BoundaryBoxUtils( const BoundaryBox<DIM> &bbox );

virtual ~BoundaryBoxUtils();

/*!
 * @brief Set boundary box.
 *
 * All utility operations refer to this box.
 */
void setBoundaryBox( const BoundaryBox<DIM> &bbox );


/*!
 * @brief Get boundary box.
 *
 * All utility operations refer to this box.
 */
const BoundaryBox<DIM> &getBoundaryBox() const;

/*!
 * @brief Get the outward direction in logical space.
 *
 * Each component of the outward direction will have
 * one of these values:
 *   - -1 if the outward direction is toward the lower indices
 *   - 0 for the direction is orthogonal to the outward direction.
 *   - 1 if the outward direction is toward the higher indices
 */
const IntVector<DIM> &getOutwardShift() const;

/*!
 * @brief Stretch box outward by the given ghost cell width.
 *
 * The number of direction affected is the same as the
 * codimension of the boundary.
 *
 * Note that the BoundaryBox is defined to be one cell wide.  The
 * output of this method is the BoundaryBox, stretched to cover the
 * given ghost cell width.  This means that if gcw is one, the output
 * is identical to the BoundaryBox.  If the gcw is zero in any
 * direction, the output will shrink to nothing in that direction.
 *
 * Return the stretched box.
 */
void stretchBoxToGhostWidth(
   Box<DIM>& box,
   const hier::IntVector<DIM> &gcw ) const;

/*!
 * @brief Extend box outward by the given amount.
 *
 * The number of directions extended is the same as the
 * codimension of the boundary.
 *
 * Return the extended box.
 */
void extendBoxOutward(
   Box<DIM> &box,
   const IntVector<DIM> &extension ) const;

/*!
 * @brief Return the normal direction.
 *
 * The normal direction is defined only for surface
 * boundaries (codimension 1).  A -1 is returned for
 * all other boundary types.
 */
int normalDir() const;

/*!
 * @brief Trim a boundary box so that it does not stick out
 * past a limiting box in direction transverse to the boundary
 * normal.
 *
 * This method affects the only box dimensions parallel to
 * the boundary.  For methods affecting other box dimensions,
 * see stretchBoxToGhostWidth().
 *
 * The boundary type of the BoundaryBox that was given to the BoundaryBoxUtils
 * constructor must be less than DIM.
 *
 * @param limit_box Box to not stick past
 *
 * @return New trimmed boundary box.
*/
hier::BoundaryBox<DIM> trimBoundaryBox(
   const hier::Box<DIM> &limit_box ) const;

/*!
 * @brief Return box describing the index space of the outer surface of
 * a boundary box.
 *
 * Define a box describing the indices of the surface of the
 * the input boundary box.  A surface is a face in 3D and an edge
 * in 2D.  These surfaces lie on the boundary itself.

 * The input boundary_box must be of type 1
 * (see hier::BoundaryBox::getBoundaryType()).

 * This is a utility function for working with the surface
 * indices coresponding to a boundary box.

 * @return a box to define the face indices corresponding to
 * boundary_box
*/
hier::Box<DIM> getSurfaceBoxFromBoundaryBox() const;

private:

/*!
 * @brief Compute the shift in the outward direction
 * (redundant data, function of boundary box) and store
 * in d_outward.
 *
 * @see getOutwardShift();
 */
void computeOutwardShift();

/*!
 * @brief Boundary box implicitly referred to by all methods.
 *
 * @see setBoundaryBox(), getBoundaryBox()
 */
BoundaryBox<DIM> d_bbox;

/*!
 * @brief Vector pointing outward from patch.
 */
IntVector<DIM> d_outward;

};


}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoundaryBoxUtils.C"
#endif
