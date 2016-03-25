//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/multiblock/MultiblockCoarsenAlgorithm.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1818 $
// Modified:	$LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
// Description:	Coarsening algorithm for data transfer between AMR levels
//
 
#ifndef included_xfer_MultiblockCoarsenAlgorithm
#define included_xfer_MultiblockCoarsenAlgorithm

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_xfer_MultiblockCoarsenSchedule
#include "MultiblockCoarsenSchedule.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_xfer_CoarsenClasses
#include "CoarsenClasses.h"
#endif

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/*!
 * @brief Class MultiblockCoarsenAlgorithm<DIM> encapsulates the AMR
 * communication pattern to coarsen data from a finer level to any coarser
 * level in a multiblock domain.  Most often, data is coarsened from the
 * interiors of source patch components on the source patch level into
 * interiors of patch components on the destination level.  If the
 * coarsening operators require ghost cells on a source component, then
 * sufficient ghost cell storage must be provided by the source patch data
 * component, and those ghost cells must be filled before calling the data
 * coarsening routines.
 *
 * Communication algorithms generally consist of three parts: an algorithm,
 * a schedule, and a patch strategy.  The algorithm describes the communication
 * between patch data items but is independent of the configuration of the
 * AMR hierarchy.  Patch data items and their associated coarsening operators
 * are registered with the algorithm.  To generate the communication
 * dependencies for a particular hierarchy configuration, the algorithm
 * generates a schedule based on the current hierarchy configuration.  This
 * schedule then performs the communication based on the registered data types
 * and their associated operators.  User-defined pre-processing and
 * post-processing are provided through the abstract patch strategy class.
 * The source patch data space is used during processing to store temporary
 * data.  Thus, the user-defined coarsening operators should operate on the
 * source space by using the patch data with those indices.
 *
 * Note that each coarsen schedule created by a coarsen algorithm remains
 * valid as long as the patches involved in the communication process do not
 * change as long as the patches involved in the communication process do not
 * change; thus, they can be used for multiple data communication cycles. 
 *
 *  Typical usage of a coarsen algorithm to perform data coarsening
 * on an AMR hierarchy involves four steps:
 *
 * -# Construct a coarsen algorithm object.
 * -# Register coarsen operations with the coarsen algorithm.  Using the
 *       registerCoarsen() methods(s), one provides source and destination
 *       patch data information, as well as spatial coarsening operators
 *       as needed.
 * -# After all operations are registered with the algorithm, one
 *       creates a communication schedule using the createSchedule()
 *       method.  This method identifies the source (fine) and destination
 *       (coarse) patch levels for data coarsening.  Note that when creating
 *       a communication schedule, a concrete instance of a
 *       MultiblockCoarsenPatchStrategy<DIM> object may be required to supply
 *       user-defined spatial data coarsening operations.
 * -# Invoke the coarsenData() method in the communication schedule to
 *       perform the data transfers.
 *
 * @see xfer::MultiblockCoarsenSchedule
 * @see xfer::CoarsenSchedule
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::CoarsenClasses
*/

template<int DIM> 
class MultiblockCoarsenAlgorithm  : public tbox::DescribedClass
{
public:
   /*!
    * Construct a coarsening algorithm and initialize its basic state.
    * Coarsening operations must be registered with this algorithm
    * before it can do anything useful.  See the registerCoarsen() routine
    * for details
    *
    * @param multiblock        pointer to the multiblock patch hierarchy
    *                          on which this coarsen algorithm with operate.
    * @param fill_coarse_data  boolean flag indicating whether coarse level data
    *                          is needed for the data coarsening operations.
    *                          By default this argument is false.  If a true
    *                          value is provided, then source data will be
    *                          filled on a temporary coarse patch level (i.e.,
    *                          copied from the actual coarse level source data)
    *                          for use in coarsening operations registered with
    *                          this algorithm.  This option should only be used
    *                          by those who specifically require this behavior
    *                          and who know how to properly process the patch
    *                          data on coarse and fine patch levels during
    *                          the coarsening process.
    */
   MultiblockCoarsenAlgorithm(
      tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
      bool fill_coarse_data = false);
 
   /*!
    * The virtual destructor for the algorithm releases all internal storage.
    */
   virtual ~MultiblockCoarsenAlgorithm<DIM>();

   /*!
    * Register a coarsening operation with the coarsening algorithm.  Data
    * from the interiors of the source component on a source (fine) patch level
    * will be coarsened into the source component of a temporary (coarse) patch
    * level and then copied into the destination component on the destination
    * (coarse) patch level.  If the coarsening operator requires data in ghost
    * cells outside of the patch interiors (i.e., a non-zero stencil width),
    * then those ghost cells must exist in the source patch data component
    * and the ghost cells must be filled with valid data on the source level
    * before a call to invoke the communication schedule.   Note that the
    * source and destination components may be the same in any case.
    *
    * Some special circumstances require that data be coarsened from the
    * ghost cell regions of a finer level and the resulting coarsened data
    * should be copied to the destination patch level.  When this is the case,
    * the optional integer vector argument should be set to the cell width, in
    * the destination (coarser) level index space, of the region around the
    * fine level where this coarsening should occur.  For example, if the
    * coarser level needs data in a region two (coarse) cells wide around the
    * boundary of the finer level, then the gcw_to_coarsen should be set to a
    * vector with all entries set to two.  Moreover, if the refinement ratio  
    * between coarse and fine levels is four in this case, then the source
    * patch data is required to have at least eight ghost cells.
    *
    * @param dst       Patch data index filled on destination level.
    * @param src       Patch data index coarsened from the source level.
    * @param opcoarsen Pointer to coarsening operator.  This may be a null
    *                  pointer.
    *
    *                  In this case, coarsening must be handled by the coarsen
    *                  patch strategy member functions.  See the comments for
    *                  MultiblockCoarsenPatchStrategy<DIM>::preprocessCoarsen() and
    *                  MultiblockCoarsenPatchStrategy<DIM>::postprocessCoarsen().
    * @param gcw_to_coarsen Integer vector ghost cell width when data should
    *                       be coarsened from ghost cell regions of the source
    *                       (finer) level into the destination (coarser) level.
    *                       By default, it is the vector of zeros indicating
    *                       that data should be coarsened from from patch
    *                       interiors on the source level.  If this argument is
    *                       used, its value should be the cell width, in the
    *                       destination (coarser) level index space, of the
    *                       region around the fine level where this coarsening
    *                       should occur.  This argument should only be
    *                       provided by those who specifically require this
    *                       special behavior and know how to properly process
    *                       the patch data on coarse and fine patch levels
    *                       during the coarsening process.
    */
   void registerCoarsen(
      const int dst,
      const int src,
      const tbox::Pointer< xfer::CoarsenOperator<DIM> > opcoarsen,
      const hier::IntVector<DIM>& gcw_to_coarsen = hier::IntVector<DIM>(0));

   /*!
    * Create a communication schedule to coarsen data from finer patch
    * levels.  This communication schedule may then be executed to perform
    * the data transfers.  This schedule creation procedure assumes that
    * the coarse level represents a ragion of coarser index space than the
    * fine level.  To avoid potentially erroneous behavior, the coarse level
    * domain should cover the domain of the fine level.
    *
    * Neither patch level can be null and when assertion checking is active,
    * passing a null level pointer will produce an unrecoverable exception.
    *
    * Note that the schedule remains valid as long as the patches on the levels
    * do not change; thus, it can be used for multiple data communication
    * cycles.
    *
    * @return Pointer to coarsen schedule that performs the data transfers.
    *
    * @param crse_level       Pointer to coarse (destination) level.
    * @param fine_level       Pointer to fine (source) level.
    * @param patch_strategy Pointer to a coarsen patch strategy that provides
    *                         user-defined coarsen operations.  If this patch
    *                         strategy is null (default state), then no
    *                         user-defined coarsen operations will be
    *                         performed.
    * @param refine_strategy  Pointer to a refine patch strategy that provides
    *                         user-defined refinement operations for specific
    *                         cases where a user-defined coarsening scheme
    *                         requires data to be filled on the coarse level
    *                         prior to the execution of the coarsening
    *                         operation.  Default is NULL. 
    */
   tbox::Pointer< MultiblockCoarsenSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > crse_level,
      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > fine_level,
      MultiblockCoarsenPatchStrategy<DIM>* patch_strategy =
            ((MultiblockCoarsenPatchStrategy<DIM>*)NULL),
      MultiblockRefinePatchStrategy<DIM>* refine_strategy =
            ((MultiblockRefinePatchStrategy<DIM>*)NULL)) const;


private:
   // These two functions are not implemented
   MultiblockCoarsenAlgorithm(const MultiblockCoarsenAlgorithm<DIM>&);
   void operator=(const MultiblockCoarsenAlgorithm<DIM>&);

   tbox::Pointer< xfer::CoarsenClasses<DIM> > d_coarsen_classes;
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > d_multiblock;
   bool d_fill_coarse_data;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockCoarsenAlgorithm.C"
#endif

#endif
