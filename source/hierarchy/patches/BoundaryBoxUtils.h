//
// File:	$URL$
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2006 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Generic utilities for boundary box calculus.
//

#ifndef included_hier_BoundaryBoxUtils
#define included_hier_BoundaryBoxUtils

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif


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
 * The number of direction extended is the same as the
 * codimension of the boundary.
 *
 * Note that the BoundaryBox is defined to be one cell
 * wide.  The output is the BoundaryBox, stretched to
 * cover the given ghost cell width.  This means that
 * if gcw is one, the output is identical to the BoundaryBox.
 * If the gcw is zero in any direction, the output will
 * shrink to nothing in that direction.
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
 * @brief Shift box inward by the given distance.
 *
 * The number of direction shifted is the same as the
 * codimension of the boundary.
 *
 * Return the shifted box.
 */
void shiftBoxInward(
   Box<DIM> &box,
   const hier::IntVector<DIM> &distance ) const;

/*!
 * @brief Return the normal direction.
 *
 * The normal direction is defined only for surface
 * boundaries (codimension 1).  A -1 is returned for
 * all other boundary types.
 */
int normalDir() const;

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

#ifndef DEBUG_NO_INLINE
// #include "BoundaryBoxUtils.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
// #include "BoundaryBoxUtils.C"
#endif
