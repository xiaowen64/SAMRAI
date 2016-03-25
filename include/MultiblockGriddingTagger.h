//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/multiblock/MultiblockGriddingTagger.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Strategy interface to user routines for refining AMR data.
//
 
#ifndef included_mesh_MultiblockGriddingTagger
#define included_mesh_MultiblockGriddingTagger

#include "SAMRAI_config.h"
#include "MultiblockRefinePatchStrategy.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace mesh {

/*!
 * @brief Class MultiblockGriddingTagger<DIM> is a concrete implementation
 * of MultiblockRefinePatchStrategy<DIMI> that is used for boundary filling
 * of patch data representing cells tagged for refinement.
 *
 * This class is needed for the calls to MultiblockRefineSchedule<DIM> in
 * the MultiblockGriddingAlgorithm<DIM>.
 *
 * This class implements the interface from MultiblockRefinePatchStrategy for
 * fillSingularityBoundaryConditions(), so that boundary conditions for
 * tag data that abuts a singularity can be properly filled.  Also
 * implemented are the interfaces for xfer::RefinePatchStrategy<DIM>, needed
 * primarily for physical boundary filling.
 *
 * @see mesh::MultiblockGriddingAlgorithm
 * @see xfer::MultiblockRefineSchedule
 * @see xfer::MultiblockRefinePatchStrategy
 * @see xfer::RefinePatchStrategy
 */

template<int DIM>
class MultiblockGriddingTagger
:
public xfer::MultiblockRefinePatchStrategy<DIM>
{
public:
   /*!
    * @brief The constructor does nothing interesting.
    */
   MultiblockGriddingTagger();
 
   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockGriddingTagger<DIM>();

   /*!
    * @brief Set the patch data index for tag data.  This routine
    * must be called with a valid cell-centered integer patch data
    * index.
    */
   virtual void setScratchTagPatchDataIndex(int buf_tag_indx);

   /*!
    * @brief Physical boundary fill
    *
    * Implementation of interface defined in xfer::RefinePatchStrategy<DIM>.
    * Fills ghost cells of patch data at physical boundaries.
    *
    * @param patch               Patch where data is stored
    * @param fill_time           Simulation time when filling occurs
    * @param ghost_width_to_fill Maximum ghost width of all data to be filled
    */
   virtual void setPhysicalBoundaryConditions(
      hier::Patch<DIM>& patch,
      const double fill_time,
      const hier::IntVector<DIM>& ghost_width_to_fill);

   /*!
    * @brief Set the ghost data at a multiblock singularity.
    *
    * Implementation of interface defined in MultiblockRefinePatchStrategy<DIM>.
    * Fills ghost cells of patch data that abut multiblock singularities.
    * The list of singularity patches contains the data from neighboring
    * blocks that also abut the singularity, and that data from the neighbors
    * is used to fill data on the local patch.
    *
    * @param patch               Local patch containing data to be filled
    * @param singularity_patches List of structures that contain data from
    *                             neighboring blocks.  See
    *                             MultiblockRefineSchedule<DIM> for more
    *                             information on struct SingularityPatch
    * @param fill_time            Simulation time when filling occurs
    * @param fill_box             All ghost data to be filled will be within
    *                             this box
    * @param boundary_box         BoundaryBox object that stores information
    *                             about the type and location of the boundary
    *                             where ghost cells will be filled
    */
   virtual void fillSingularityBoundaryConditions(
      hier::Patch<DIM>& patch,
      tbox::List<typename xfer::MultiblockRefineSchedule<DIM>::SingularityPatch>&
         singularity_patches, 
      const double fill_time,
      const hier::Box<DIM>& fill_box,
      const hier::BoundaryBox<DIM>& boundary_box);

   /*!
    * @brief Return maximum stencil width needed for user-defined
    * data interpolation operations.  This is needed to
    * determine the correct interpolation data dependencies.
    *
    * Always returns an IntVector of ones, because that is the maximum
    * stencil needed for the operations in MultiblockGriddingAlgorithm
    */
   virtual hier::IntVector<DIM> getRefineOpStencilWidth() const
   {
      return (hier::IntVector<DIM>(1));
   }

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
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
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
      const hier::IntVector<DIM>& ratio);

private:
   /* 
    * Patch data index for 
    */
   int d_buf_tag_indx;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockGriddingTagger.C"
#endif

#endif
