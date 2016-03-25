//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/multiblock/MultiblockRefinePatchStrategy.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Strategy interface to user routines for refining AMR data.
//
 
#ifndef included_xfer_MultiblockRefinePatchStrategy
#define included_xfer_MultiblockRefinePatchStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_xfer_MultiblockRefineSchedule
#include "MultiblockRefineSchedule.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_xfer_CoarsenPatchStrategy
#include "CoarsenPatchStrategy.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif

#ifndef MULTIBLOCK_UNDEFINED_BLOCK_NUMBER
#define MULTIBLOCK_UNDEFINED_BLOCK_NUMBER (-1)
#endif

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class MultiblockRefinePatchStrategy is a virtual base class
 * that provides interfaces for users to problem-specific routines related
 * to issues that arise in multiblock domains, particularly the filling
 * of boundary conditions around a singularity.
 */

template<int DIM>
class MultiblockRefinePatchStrategy
:
public virtual tbox::DescribedClass,
public virtual xfer::RefinePatchStrategy<DIM>
{
public:
   /*!
    * The constructor for patch strategy does nothing interesting.
    */
   MultiblockRefinePatchStrategy();
 
   /*!
    * The virtual destructor for refine strategy does nothing interesting.
    */
   virtual ~MultiblockRefinePatchStrategy<DIM>();

   /*!
    * @brief Set the physical boundary conditions.
    *
    * @param patch The patch containing the data to be filled
    * @param fill_time Simulation time at which data is filled
    * @param ghost_width_to_fill maximum number of ghost cells to fill
    */
   virtual void setPhysicalBoundaryConditions(
      hier::Patch<DIM>& patch,
      const double fill_time,
      const hier::IntVector<DIM>& ghost_width_to_fill);

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
   virtual void fillSingularityBoundaryConditions(
      hier::Patch<DIM>& patch,
      tbox::List<typename MultiblockRefineSchedule<DIM>::SingularityPatch>&
         singularity_patches, 
      const double fill_time,
      const hier::Box<DIM>& fill_box,
      const hier::BoundaryBox<DIM>& boundary_box) = 0;

   /*!
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  This is needed to
    * determine the correct interpolation data dependencies.
    */
   virtual hier::IntVector<DIM> getRefineOpStencilWidth() const = 0;

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before standard refining operations (expressed using
    * concrete subclasses of the xfer::RefineOperator<DIM> base class).
    * The preprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm<DIM> class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
      return;
   }

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before standard refining operations (expressed using
    * concrete subclasses of the xfer::RefineOperator<DIM> base class).
    * The postprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm<DIM> class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void postprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
      return;
   }

   /*
    * During the fillData operation by the MultiblockRefineSchedule, there are 
    * times where it is necssary to fill a temporary coarsened "scratch"
    * level.  This method allows the schedule to set the right conditions
    * in this strategy class.
    */
   virtual void setFillingCoarseScratch(const bool filling_coarse_scratch)
   {
      d_filling_coarse_scratch = filling_coarse_scratch;
   }

   /*
    * Set the multiblock block number.
    */
   virtual void setBlockNumber(const int block_number)
   {
      d_block_number = block_number;
   }

   /*
    * Get the multiblock block number.
    */
   virtual int getBlockNumber()
   {
      return(d_block_number);
   }

   /*
    * Clear the multiblock block number.
    */
   virtual void clearBlockNumber()
   {
      d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
   }

protected:

   bool d_filling_coarse_scratch;
   int  d_block_number;

private:

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockRefinePatchStrategy.C"
#endif

#endif
