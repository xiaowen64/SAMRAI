//
// File:	RefinePatchStrategy.h
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Strategy interface to user routines for refining AMR data.
//
 
#ifndef included_xfer_RefinePatchStrategy
#define included_xfer_RefinePatchStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class RefinePatchStrategy<DIM> is an abstract base class that
 * defines a Strategy pattern interface for operations that a user must 
 * supply for application-specific refining of data between levels in an
 * AMR patch hierarchy and filling of physical boundary values.   
 * A concrete subclass must define four member functions to perform 
 * the following tasks:
 *
 * - define maximum stencil width for user-defined refine operations
 * - fill physical boundary conditions
 * - preprocess the interpolation
 * - postprocess the interpolation
 *
 * Note that the preprocess member function is called before standard data refine
 * using RefineOperators and the postprocessor member function is called afterwards.
 * 
 * There are two versions of the preprocess and postprocess functions.  The
 * default abstract function only takes a single box.  The user may also
 * over-ride preprocess and postprocess functions that take a box list and
 * process an entire patch at one time.  By default, the box list version
 * loops over all of the boxes in the box list and calls the other version.
 *
 * @see xfer::RefineAlgorithm
 * @see xfer::RefineSchedule
 */
 
template<int DIM> class RefinePatchStrategy
:
public virtual tbox::DescribedClass
{
public:
   /*!
    * The constructor for patch strategy does nothing interesting.
    */
   RefinePatchStrategy();
 
   /*!
    * The virtual destructor for refine strategy does nothing interesting.
    */
   virtual ~RefinePatchStrategy();

   /*!
    * Pure virtual function to set data associated with the given list of patch data
    * indices at patch boundaries that intersect the physical domain boundary.  The
    * specific boundary conditions are determined by the user.  The patch data
    * components set in this routine correspond to the "scratch" components specified
    * in calls to the registerRefine() function in the RefineAlgorithm<DIM> class.
    *
    * @param patch                hier::Patch on which to fill boundary data.
    * @param fill_time            Double simulation time for boundary filling.
    * @param ghost_width_to_fill  Integer vector describing maximum ghost width
    *                             to fill over all registered scratch components.
    */
   virtual void setPhysicalBoundaryConditions(
      hier::Patch<DIM>& patch,
      const double fill_time,
      const hier::IntVector<DIM>& ghost_width_to_fill) = 0;

   /*!
    * Pure virtual function to return maximum stencil width needed over
    * user-defined data interpolation operations.  This is needed to
    * determine the correct interpolation data dependencies.
    */
   virtual hier::IntVector<DIM> getRefineOpStencilWidth() const = 0; 

   /*!
    * Pure virtual function to perform user-defined preprocess data refine 
    * operations.  This member function is called before standard refine 
    * operations (expressed using concrete subclasses of the RefineOperator<DIM> 
    * base class). The preprocess function must refine data from the scratch 
    * components of the coarse patch into the scratch components of the fine 
    * patch on the specified fine box region.  Recall that the scratch components 
    * are specified in calls to the registerRefine() function in the 
    * RefineAlgorithm<DIM> class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    hier::Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio) = 0;

   /*!
    * Pure virtual function to perform user-defined preprocess data refine
    * operations.  This member function is called after standard refine
    * operations (expressed using concrete subclasses of the RefineOperator<DIM>
    * base class). The postprocess function must refine data from the scratch 
    * components of the coarse patch into the scratch components of the fine patch 
    * on the specified fine box region.  Recall that the scratch components are 
    * specified in calls to the registerRefine() function in the 
    * RefineAlgorithm<DIM> class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    hier::Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void postprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio) = 0;

   /*!
    * Virtual function to perform user-defined refine operations.  This member
    * function is called before standard refining operations (expressed using
    * concrete subclasses of the RefineOperator<DIM> base class).  The preprocess
    * function must refine data from the scratch components of the coarse patch
    * into the scratch components of the fine patch on the specified fine box
    * regions.
    *
    * Typically, only the pure virtual members of this class are implemented in
    * user-defined subclasses of this base class.   This version of the preprocess
    * function operates on an entire box list.  By default, this version simply
    * loops over the box list and calls the single-box version, which is a pure
    * virtual function.
    * 
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_boxes  tbox::List of box regions on fine patch into which data 
    *                    is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessRefineBoxes(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::BoxList<DIM>& fine_boxes,
      const hier::IntVector<DIM>& ratio);

   /*!
    * Virtual function to perform user-defined refine operations.  This member
    * function is called after standard refining operations (expressed using
    * concrete subclasses of the RefineOperator<DIM> base class).  The postprocess
    * function must refine data from the scratch components of the coarse patch
    * into the scratch components of the fine patch on the specified fine box
    * regions.
    *
    * Typically, only the pure virtual members of this class are implemented in
    * user-defined subclasses of this base class.   This version of the postprocess
    * function operates on an entire box list.  By default, this version simply
    * loops over the box list and calls the single-box version, which is a pure
    * virtual function.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_boxes  tbox::List of box regions on fine patch into which data
    *                    is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void postprocessRefineBoxes(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::BoxList<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefinePatchStrategy.C"
#endif
