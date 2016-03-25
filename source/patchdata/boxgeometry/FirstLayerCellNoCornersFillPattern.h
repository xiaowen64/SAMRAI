//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	First layer cell fill pattern class
//
 
#ifndef included_pdat_FirstLayerCellNoCornersFillPattern
#define included_pdat_FirstLayerCellNoCornersFillPattern

#include "SAMRAI_config.h"

#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "VariableFillPattern.h"

namespace SAMRAI {
    namespace pdat {


/*!
 * Class FirstLayerCellNoCornersFillPattern is a concrete implementation of
 * the abstract base class VariableFillPattern.  It is used to calculate
 * overlaps according to a pattern which limits overlaps to the cell-centered
 * ghost region of width 1 surrounding a patch, excluding all edges and
 * corners.
 */

template<int DIM>
class FirstLayerCellNoCornersFillPattern :
public xfer::VariableFillPattern<DIM>
{
public:

   /*!
    * Constructor
    *
    * @param dim     Dimension
    */
   FirstLayerCellNoCornersFillPattern();

   /*!
    * Destructor
    */
   virtual ~FirstLayerCellNoCornersFillPattern();

   /*!
    * Calculate overlaps between the destination and source geometries
    * according to the desired pattern.  This will return the portion
    * of the intersection of the geometries that lies in the ghost region
    * of width 1 surrounding the patch, excluding all edges and corners.
    * The patch is identified by the argument dst_patch_box.
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
    * Returns the stencil width of 1 in all directions.
    */
   hier::IntVector<DIM>& getStencilWidth();

   /*!
    * Returns a string name identifier
    * "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN".
    */
   const std::string& getPatternName() const;

private:
   FirstLayerCellNoCornersFillPattern(
      const FirstLayerCellNoCornersFillPattern<DIM>&);    // not implemented
   void operator=(const FirstLayerCellNoCornersFillPattern<DIM>&);         // not implemented

   hier::IntVector<DIM> d_stencil_width;
   static std::string s_name_id;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FirstLayerCellNoCornersFillPattern.C"
#endif

