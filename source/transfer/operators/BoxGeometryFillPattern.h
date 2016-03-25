//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Default fill pattern class
//
 
#ifndef included_xfer_BoxGeometryFillPattern
#define included_xfer_BoxGeometryFillPattern

#include "SAMRAI_config.h"

#include "VariableFillPattern.h"

namespace SAMRAI {
    namespace xfer {


/*!
 * Class BoxGeometryFillPattern is a default implementation of
 * the abstract base class VariableFillPattern.  It is used to calculate
 * overlaps that consist of the intersections between boxes, including
 * the full ghost regions.
 */
template<int DIM>
class BoxGeometryFillPattern :
public VariableFillPattern<DIM>
{
public:

   /*!
    * Default constructor
    */
   BoxGeometryFillPattern();

   /*!
    * Destructor
    */
   virtual ~BoxGeometryFillPattern();

   /*!
    * Calculate overlaps between the destination and source geometries
    * using the geometries' own overlap calculation methods.
    *
    * @param dst_geometry    geometry object for destination box
    * @param src_geometry    geometry object for source box
    * @param dst_patch_box   box for the destination patch
    * @param src_mask        the source mask, the box resulting from shifting
    *                        the source box
    * @param overwrite_interior  controls whether or not to include the
    *                            destination box interior in the overlap
    * @param src_offset      the offset between source and destination
    *                        index space.  src + src_offset = dst
    *
    * @return                Pointer to the calculated overlap object
    */
   tbox::Pointer< hier::BoxOverlap<DIM> >
   calculateOverlap(const hier::BoxGeometry<DIM>& dst_geometry,
                    const hier::BoxGeometry<DIM>& src_geometry,
                    const hier::Box<DIM>& dst_patch_box,
                    const hier::Box<DIM>& src_mask,
                    const bool overwrite_interior,
                    const hier::IntVector<DIM>& src_offset) const;

   /*!
    * Implementation of interface to get stencil width of a
    * VariableFillPattern.  For this class BoxGeometryFillPattern, this
    * method should never be called, since overlaps are computed based
    * on box geometry objects and not on any stencil.  An error will result
    * if this method is invoked.
    */
   hier::IntVector<DIM>& getStencilWidth();

   /*!
    * Returns a string name identifier "BOX_GEOMETRY_FILL_PATTERN".
    */
   const std::string& getPatternName() const;

private:
   BoxGeometryFillPattern(
      const BoxGeometryFillPattern<DIM>&);    // not implemented
   void operator=(const BoxGeometryFillPattern<DIM>&);         // not implemented

   static std::string s_name_id;
   static hier::IntVector<DIM> s_stencil_width;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGeometryFillPattern.C"
#endif

