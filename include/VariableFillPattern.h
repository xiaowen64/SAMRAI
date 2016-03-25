//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Abstract fill pattern class to provide interface for stencils
//
 
#ifndef included_xfer_VariableFillPattern
#define included_xfer_VariableFillPattern

#include "SAMRAI_config.h"

#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <string>

namespace SAMRAI {
    namespace xfer {


/*!
 * Class VariableFillPattern is an abstract base class that provides an 
 * interface to create objects that can calculate overlaps which correspond
 * to a specific stencil.
 */
template<int DIM>
class VariableFillPattern : public tbox::DescribedClass
{
public:

   /*!
    * Default constructor 
    */
   VariableFillPattern();

   /*!
    * Destructor
    */
   virtual ~VariableFillPattern();

   /*!
    * This pure virtual method provides an interface to calculate overlaps
    * between the destination and source geometries.
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
   virtual tbox::Pointer< hier::BoxOverlap<DIM> >
   calculateOverlap(const hier::BoxGeometry<DIM>& dst_geometry,
                    const hier::BoxGeometry<DIM>& src_geometry,
                    const hier::Box<DIM>& dst_patch_box,
                    const hier::Box<DIM>& src_mask,
                    const bool overwrite_interior,
                    const hier::IntVector<DIM>& src_offset) const = 0;

   /*!
    * Return the maximum ghost width of the boundary stencil.  The default
    * implementation throws an error.
    */
   virtual hier::IntVector<DIM>& getStencilWidth() = 0;

   /*!
    * Return a string name identifying the concrete subclass.
    */
   virtual const std::string& getPatternName() const = 0;

private:
   VariableFillPattern(const VariableFillPattern<DIM>&);    // not implemented
   void operator=(const VariableFillPattern<DIM>&);         // not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "VariableFillPattern.C"
#endif
