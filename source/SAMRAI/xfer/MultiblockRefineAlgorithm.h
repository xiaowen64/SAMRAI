/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   class to manage multiblocks 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefineAlgorithm
#define included_xfer_MultiblockRefineAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/MultiblockRefineSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/MultiblockRefinePatchStrategy.h"
#include "SAMRAI/tbox/Pointer.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class MultiblockRefineAlgorithm is an extension of the
 * concept of xfer::RefineAlgorithm to be used in applications that require
 * a multiblock domain.
 *
 * This class requires a pointer to a previously-constructed
 * xfer::RefineAlgorithm.  It contains four routines that create
 * a MultiblockRefineSchedule object, each analagous to the
 * createSchedule() routines in xfer::RefineAlgorithm.
 *
 * @see PatchHierarchy
 * @see MultiblockRefineSchedule
 * @see MultiblockRefineAlgorithm
 * @see RefineAlgorithm
 */

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
    * @param hierarchy     The Multiblock patch hierarchy on which this
    *                      algorithm will operate.
    * @param dim           Dimension for the object
    */
   explicit MultiblockRefineAlgorithm(
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   ~MultiblockRefineAlgorithm();

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
   void
   registerRefine(
      const int dst,
      const int src,
      const int scratch,
      tbox::Pointer<hier::RefineOperator> oprefine,
      tbox::Pointer<VariableFillPattern> var_fill_pattern =
         (tbox::Pointer<BoxGeometryVariableFillPattern>)NULL);

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
    * in an unrecoverable assertion.
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
   void
   registerRefine(
      const int dst,
      const int src,
      const int src_told,
      const int src_tnew,
      const int scratch,
      tbox::Pointer<hier::RefineOperator> oprefine,
      tbox::Pointer<hier::TimeInterpolateOperator> optime,
      tbox::Pointer<VariableFillPattern> var_fill_pattern =
         (tbox::Pointer<BoxGeometryVariableFillPattern>)NULL);

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
    * @param level           Pointer to PatchLevel which contains
    *                        PatchLevels on which interpatch transfers
    *                        occur.  This pointer cannot be null.
    * @param patch_strategy   Pointer to a multiblock patch strategy that
    *                        provides user-defined singularity boundary filling
    *                        operations.
    * @param transaction_factory  Refine transaction factory.
    */
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<hier::PatchLevel> level,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

   /*!
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<PatchLevelFillPattern> fill_pattern,
      tbox::Pointer<hier::PatchLevel> level,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

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
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

   /*!
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<PatchLevelFillPattern> fill_pattern,
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

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
    * @param level           Pointer to destination PatchLevel.  This
    *                        This pointer cannot be null.
    * @param next_coarser_level Integer number of next coarser
    *                           PatchLevel relative to the
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
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<hier::PatchLevel> level,
      const int next_coarser_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

   /*!
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<PatchLevelFillPattern> fill_pattern,
      tbox::Pointer<hier::PatchLevel> level,
      const int next_coarser_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;
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
    *                           PatchLevel relative to the
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
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      const int next_coarser_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

   /*!
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer<MultiblockRefineSchedule>
   createSchedule(
      tbox::Pointer<PatchLevelFillPattern> fill_pattern,
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      const int next_coarser_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      MultiblockRefinePatchStrategy* patch_strategy =
         ((MultiblockRefinePatchStrategy *)NULL),
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory =
         tbox::Pointer<xfer::RefineTransactionFactory>(NULL)) const;

   /*!
    * @brief Given a previously-generated refine schedule, reconfigure it to
    * peform the communication operations registered with this refine algorithm
    * object.  That is, the schedule will be transformed so that it will
    * funcions as though this refine algorithm created it.  Note that the set
    * of operations registered with this refine algorithm must be essentially
    * the same as those registered with the refine algorithm that created the
    * schedule originallyl.  That is, the number of operations registered must
    * be the same and the source, destination, scratch patch data items and
    * operators for each operation must have identical characteristics (i.e.,
    * data centering, ghost cell widths, stencil requirements, etc.).
    * However, the specific source, destination, scratch patch data ids and
    * operators can be different.  Detailed and fairly complete error checking
    * is performed when this routine is called to prevent potential errors or
    * unexpected behavior.
    *
    * @param schedule  tbox::Pointer to refine schedule, which cannot be null.
    */
   void
   resetSchedule(
      tbox::Pointer<xfer::MultiblockRefineSchedule> schedule) const;

   /*!
    * @brief Disable the use of singularity patch functionality in all
    * schedules created for this algorithm.
    *
    * This should only be used in specific cases when it is desired that
    * no ghost data be communicated or filled at enhanced or reduced
    * connectivity patch boundaries.
    */
   void
   disableSingularityPatches()
   {
      d_enable_singularity_patches = false;
   }

private:
   tbox::Pointer<xfer::RefineAlgorithm> d_single_block_refine_alg;
   tbox::Pointer<hier::PatchHierarchy> d_multiblock_hierarchy;

   bool d_enable_singularity_patches;
};

}
}

#endif
