//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/multiblock/MultiblockRefineAlgorithm.h $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 3153 $
// Modified:    $LastChangedDate: 2009-04-21 17:12:47 -0700 (Tue, 21 Apr 2009) $
// Description: class to manage multiblocks
//

#ifndef included_xfer_MultiblockRefineAlgorithm
#define included_xfer_MultiblockRefineAlgorithm

#include "SAMRAI_config.h"

#include "MultiblockPatchHierarchy.h"
#include "MultiblockPatchLevel.h"
#include "MultiblockRefineSchedule.h"
#include "tbox/Pointer.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class MultiblockRefineAlgorithm<DIM> is an extension of the
 * concept of xfer::RefineAlgorithm<DIM> to be used in applications that require
 * a multiblock domain.
 *
 * This class requires a pointer to a previously-constructed
 * xfer::RefineAlgorithm<DIM>.  It contains four routines that create
 * a MultiblockRefineSchedule<DIM> object, each analagous to the
 * createSchedule() routines in xfer::RefineAlgorithm<DIM>.
 * 
 * @see MultiblockPatchHierarchy
 * @see MultiblockRefineSchedule
 * @see MultiblockRefineAlgorithm
 * @see xfer::RefineAlgorithm
 */ 

template<int DIM>
class MultiblockRefineAlgorithm
{

public:

   /*!
    * @brief MultiblockRefineAlgorithm constructor
    *
    * Constructor for MultiblockRefineAlgorithm takes a pointer to
    * an xfer::RefineAlgorithm, that is used to create communication
    * schedules.  Refinement operations should be registered with the
    * RefineAlgorithm before this class will do anything useful.
    *
    * @param refine_alg    Pointer to a refine algorithm that should
    *                      already be constructed.
    * @param multiblock    The Multiblock patch hierarchy on which this
    *                      algorithm will operate.
    */
   MultiblockRefineAlgorithm(
      tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock);

   /*!
    * @brief Destructor
    */
   ~MultiblockRefineAlgorithm<DIM>();

   /*!
    * Register a refine operation with the refine algorithm object.  This
    * routine does not support time interpolation.  Data values will be moved
    * from the source component to the destination component using scratch
    * component as a temporary work space.  The scratch component must have
    * sufficient ghost cells to cover the required operator stencil width and
    * any needed physical boundary ghost cells.
    *
    * @param dst       Patch data index filled on the destination level.
    * @param src       Source patch data index on the source level.
    * @param scratch   Patch data index used as a temporary work space.
    * @param oprefine  Pointer to refinement operator.  This may be a null
    *                  pointer.
    *                  In this case, refinement must be handled by the refine
    *                  patch strategy member functions.  See the comments for
    *                  xfer::RefinePatchStrategy<DIM>::preprocessRefine() and
    *                  xfer::RefinePatchStrategy<DIM>::postprocessRefine().
    */
   void registerRefine(
      const int dst,
      const int src,
      const int scratch,
      tbox::Pointer< xfer::RefineOperator<DIM> > oprefine);

   /*!
    * Register a refine operation with the refine algorithm object.  This
    * routine supports time interpolation.  Time interpolation will take place
    * between the old and new source data components on coarser levels.  On
    * the destination level, data will be moved from the source component to
    * the destination component using scratch component as a temporary work
    * space.  The scratch component must have sufficient ghost cells to cover
    * the required operator stencil width and any needed physical boundary
    * ghost cells.  The time interpolation operator cannot be null.  When
    * assertion checking is active,  passing in a null pointer will result
    * in an unrecoverable exception.
    *
    * @param dst       Patch data index filled on the destination level.
    * @param src       Source patch data index on source level and for the
    *                  patch interiors on the destination level.
    * @param src_told  Old source patch data index for time interpolation.
    * @param src_tnew  New source patch data index for time interpolation.
    * @param scratch   Patch data index used as a temporary work space.
    * @param oprefine  Pointer to refinement operator.  This may be a null
    *                  pointer.  In this case, refinement must be handled by
    *                  the refine patch strategy member functions.  See the
    *                  comments for
    *                  xfer::RefinePatchStrategy<DIM>::preprocessRefine() and
    *                  xfer::RefinePatchStrategy<DIM>::postprocessRefine().
    * @param optime    Pointer to time interpolation operator.  This pointer
    *                  may not be null.
    */
   void registerRefine(
      const int dst,
      const int src,
      const int src_told,
      const int src_tnew,
      const int scratch,
      tbox::Pointer< xfer::RefineOperator<DIM> > oprefine,
      tbox::Pointer< xfer::TimeInterpolateOperator<DIM> > optime);

   /*!
    * @brief Create a communication schedule that copies data from the
    * interiors of the source components into the interior and boundary
    * cells of the destination components on the same level where those
    * sources and destinations overlap.
    *
    * Neither time nor spatial interpolation is performed.
    *
    * @return Pointer to MultiblockRefineSchedule that performs the data
    *         transfers.
    *
    * @param level           Pointer to MultiblockPatchLevel which contains
    *                        PatchLevels on which interpatch transfers
    *                        occur.  This pointer cannot be null.
    * @param patch_strategy   Pointer to a multiblock patch strategy that
    *                        provides user-defined singularity boundary filling
    *                        operations.
    * @param transaction_factory  Refine transaction factory.
    */
   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

   /*!
    * @brief Create a communication schedule that copies data from the
    * interiors of the source components on a source level into the interior
    * and boundary cells of the destination components on a destination level
    * where those sources and destinations overlap.
    *
    * Neither time nor spatial interpolation is performed.
    *
    * @return Pointer to MultiblockRefineSchedule that performs the data
    *         transfers.
    *
    * @param dst_level       Pointer to destination level; cannot be null.
    * @param src_level       Pointer to source level; cannot be null.
    * @param patch_strategy   Pointer to a multiblock patch strategy that
    *                        provides user-defined singularity boundary filling
    *                        operations.
    * @param transaction_factory  Refine transaction factory.
    */
   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

   /*!
    * @brief Create a communication schedule that copies data from the
    * interiors of the source components into the interior and boundary
    * cells of the destination components on the same level where those
    * sources and destinations overlap.  Data is moved from old and new
    * sources on coarser levels when time interpolation is needed and from
    * the source components on the patch level into the destination
    * components.
    *
    * @return Pointer to MultiblockRefineSchedule that performs the data
    *         transfers.
    *
    * @param level           Pointer to destination MultiblockPatchLevel.  This
    *                        This pointer cannot be null.
    * @param next_coarser_level Integer number of next coarser
    *                           MultiblockPatchLevel relative to the
    *                           destination level.  Note that when the
    *                           destination level has number zero (i.e., the
    *                           coarsest level), this value should be < 0.
    * @param multiblock      Pointer to Multiblock object which contains
    *                        all of the hierarchies from which data to fill
    *                        level should come.
    * @param patch_strategy   Pointer to a multiblock patch strategy that
    *                        provides user-defined singularity boundary filling
    *                        operations.
    * @param transaction_factory  Refine transaction factory.
    */
   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
      const int next_coarser_level,
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
      const int next_coarser_level,
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

   /*!
    * @brief Create a communication schedule that moves data from the
    * interiors of the source level and coarser hierarchy levels into the
    * interior and boundary cells of the destination components on the
    * destination level where those sources and destinations overlap.  Data is
    * moved from old and new sources on coarser levels when time interpolation
    * is needed and from the source components on the source level into the
    * destination components.
    *
    * @return Pointer to MultiblockRefineSchedule that performs the data
    *         transfers.
    * @param dst_level       Pointer to destination level; cannot be null.
    * @param src_level       Pointer to source level; cannot be null.
    * @param next_coarser_level Integer number of next coarser
    *                           MultiblockPatchLevel relative to the
    *                           destination level.  Note that when the
    *                           destination level has number zero (i.e., the
    *                           coarsest level), this value should be < 0.
    * @param multiblock      Pointer to Multiblock object which contains
    *                        all of the hierarchies where the PatchLevels
    *                        exist.
    * @param patch_strategy   Pointer to a multiblock patch strategy that
    *                        provides user-defined singularity boundary filling
    *                        operations.
    * @param transaction_factory  Refine transaction factory.
    */ 
   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
      const int next_coarser_level,
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

   tbox::Pointer< MultiblockRefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
      const int next_coarser_level,
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
      MultiblockRefinePatchStrategy<DIM>* patch_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL) const;

private:

   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_single_block_refine_alg;
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > d_multiblock_hierarchy;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockRefineAlgorithm.C"
#endif

#endif
