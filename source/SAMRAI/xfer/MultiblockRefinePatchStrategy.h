/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data. 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefinePatchStrategy
#define included_xfer_MultiblockRefinePatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class MultiblockRefinePatchStrategy is a virtual base class
 * that provides interfaces for users to problem-specific routines related
 * to issues that arise in multiblock domains, particularly the filling
 * of boundary conditions around a singularity.
 */

class MultiblockRefinePatchStrategy:
   public virtual tbox::DescribedClass,
   public xfer::RefinePatchStrategy
{
public:
   /*!
    * The constructor for patch strategy does nothing interesting.
    */
   explicit MultiblockRefinePatchStrategy(
      const tbox::Dimension& dim);

   /*!
    * The virtual destructor for refine strategy does nothing interesting.
    */
   virtual ~MultiblockRefinePatchStrategy();

   /*!
    * @brief Set the physical boundary conditions.
    *
    * @param patch The patch containing the data to be filled
    * @param fill_time Simulation time at which data is filled
    * @param ghost_width_to_fill maximum number of ghost cells to fill
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);

   /*!
    * @brief Set the ghost data at a multiblock singularity.
    *
    * @param patch The patch containing the data to be filled
    * @param singularity_patches structures from refine schedule that contain
    *                            patches that hold data from neighboring
    *                            blocks around the singularity
    * @param fill_time Simulation time at which data is filled
    * @param fill_box Box covering maximum amount of ghost cells to be filled
    * @param boundary_box BoundaryBox describing location of singularity in
    *                     relation to patch
    */
   virtual void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      tbox::List<tbox::Pointer<hier::Patch> >&
      singularity_patches,
      const double fill_time,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box) = 0;

   /*!
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  This is needed to
    * determine the correct interpolation data dependencies.
    */
   virtual hier::IntVector
   getRefineOpStencilWidth() const = 0;

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before standard refining operations (expressed using
    * concrete subclasses of the hier::RefineOperator base class).
    * The preprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /*!
    * Perform user-defined refining operations.  This member function
    * is called after standard refining operations (expressed using
    * concrete subclasses of the hier::RefineOperator base class).
    * The postprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /*
    * During the fillData operation by the MultiblockRefineSchedule, there are
    * times where it is necssary to fill a temporary coarsened "scratch"
    * level.  This method allows the schedule to set the right conditions
    * in this strategy class.
    */
   virtual void setFillingCoarseScratch(
      const bool filling_coarse_scratch)
   {
      d_filling_coarse_scratch = filling_coarse_scratch;
   }

protected:
   bool d_filling_coarse_scratch;

private:

};

}
}

#endif
