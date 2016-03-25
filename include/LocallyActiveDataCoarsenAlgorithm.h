//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenAlgorithm.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Coarsening algorithm for locally-active data transfer between AMR levels
//
 
#ifndef included_xfer_LocallyActiveDataCoarsenAlgorithm
#define included_xfer_LocallyActiveDataCoarsenAlgorithm

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "LocallyActiveDataPatchLevelManager.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "CoarsenClasses.h"
#include "CoarsenOperator.h"
#include "LocallyActiveDataCoarsenSchedule.h"
#include "LocallyActiveDataCoarsenPatchStrategy.h"
#include "LocallyActiveDataCoarsenTransactionFactory.h"

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/*!
 * @brief Class LocallyActiveDataCoarsenAlgorithm<DIM> encapsulates the 
 * AMR communication pattern to coarsen locally-active data from a finer 
 * level to a coarser level where the data exists on both the coarse and 
 * fine levels.  Most often, data is coarsened from the interiors of 
 * source patch components on the source patch level into interiors of 
 * destination patch components on the destination level.  See comments 
 * for the coarsen algorithm constructor for variations that are possible 
 * for (adventurous?) users. If the coarsening operators require ghost cells 
 * on a source component, then sufficient ghost cell storage must be provided 
 * by the source patch data component, and those ghost cells must be filled 
 * before calling the data coarsening routines.
 *
 * Note that this algorithm class is similar to the CoarsenAlgorithm<DIM>
 * class.  However, unlike the standard coarsen algorithm class, this one
 * does not support resetting communication schedules.
 *
 * Communication algorithms generally consist of three parts: an algorithm,
 * a schedule, and a patch strategy.  The algorithm describes the communication
 * between patch data items but is independent of the configuration of the
 * AMR hierarchy.  Patch data items and their associated coarsening operators
 * are registered with the algorithm.  To generate the communication
 * dependencies for a particular hierarchy configuration, the algorithm
 * generates a schedule based on the current hierarchy configuration.  This 
 * schedule then performs the communication based on the registered data types 
 * and their associated operators.  User-defined pre-processing and post-processing
 * and provided through the abstract patch strategy class.
 *
 * The source patch data space is used during processing to store temporary
 * data.  Thus, the user-defined coarsening operators should operate on the
 * source space by using the patch data with those indices.
 *
 * Note that each coarsen schedule created by a coarsen algorithm remains valid as
 * long as the patches involved in the communication process do not change; thus,
 * they can be used for multiple data communication cycles.
 *
 * Typical usage of a coarsen algorithm to perform data coarsening 
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
 *       LocallyActiveDataCoarsenPatchStrategy<DIM> object may be required to supply 
 *       user-defined spatial data coarsening operations.
 * -# Invoke the coarsenData() method in the communication schedule to
 *       perform the data transfers.
 *
 * @see xfer::LocallyActiveDataCoarsenSchedule
 * @see xfer::LocallyActiveDataCoarsenPatchStrategy
 * @see xfer::CoarsenClasses
 */

template<int DIM> 
class LocallyActiveDataCoarsenAlgorithm : public tbox::DescribedClass
{
public:
   /*!
    * Construct a coarsening algorithm and initialize its basic state.  
    * Coarsening operations must be registered with this algorithm 
    * before it can do anything useful.  See the registerCoarsen() routine
    * for details
    *
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
    *
    * Note that currently the described behavior above for the case that
    * fill_coarse_data = true is disabled.
    */
   LocallyActiveDataCoarsenAlgorithm(bool fill_coarse_data = false);
 
   /*!
    * The virtual destructor for the algorithm releases all internal storage.
    */
   virtual ~LocallyActiveDataCoarsenAlgorithm<DIM>();

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
    * fine level where this coarsening should occur.  For example, if the coarser
    * level needs data in a region two (coarse) cells wide around the boundary 
    * of the finer level, then the gcw_to_coarsen should be set to a vector with
    * all entries set to two.  Moreover, if the refinement ratio between coarse
    * and fine levels is four in this case, then the source patch data is required
    * to have at least eight ghost cells.
    * 
    * @param dst       Patch data index filled on destination level.
    * @param src       Patch data index coarsened from the source level.
    * @param opcoarsen Pointer to coarsening operator.  This may be a null pointer.
    *                  In this case, coarsening must be handled by the coarsen patch 
    *                  strategy member functions.  See the comments for 
    *                  preprocessCoarsen() and postprocessCoarsen() functions
    *                  in the LocallyActiveDataCoarsenPatchStrategy<DIM> class.
    * @param gcw_to_coarsen Integer vector ghost cell width when data should 
    *                       be coarsened from ghost cell regions of the source 
    *                       (finer) level into the destination (coarser) level.  
    *                       By default, it is the vector of zeros indicating that 
    *                       data should be coarsened from from patch interiors on
    *                       the source level.  If this argument is used, its value 
    *                       should be the cell width, in the destination (coarser) 
    *                       level index space, of the region around the fine level 
    *                       where this coarsening should occur.  This argument 
    *                       should only be provided by those who specifically require 
    *                       this special behavior and know how to properly process 
    *                       the patch data on coarse and fine patch levels during the 
    *                       coarsening process.
    */
   void registerCoarsen(
      const int dst,
      const int src,
      const tbox::Pointer< xfer::CoarsenOperator<DIM> > opcoarsen,
      const hier::IntVector<DIM>& gcw_to_coarsen = hier::IntVector<DIM>(0));

   /*!
    * Create a communication schedule to coarsen data from the  patch 
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
    * do not change; thus, it can be used for multiple data communication cycles.
    *
    * @return Pointer to coarsen schedule that performs the data transfers. 
    * 
    * @param crse_level     Pointer to coarse (destination) level; cannot be null.
    * @param crse_level_mgr Pointer to coarse level data manager; cannot be null.
    * @param fine_level     Pointer to fine (source) level; cannot be null.
    * @param fine_level_mgr Pointer to fine level data manager; cannot be null.
    * @param patch_strategy Pointer to a coarsen patch strategy that provides 
    *                       user-defined coarsen operations.  If this patch 
    *                       strategy is null (default state), then no 
    *                       user-defined coarsen operations will be performed.  
    * @param transaction_factory Optional tbox::Pointer to a coarsen transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardLocallyActiveDataCoarsenTransactionFactory 
    *                            object will be used.
    */
   tbox::Pointer< xfer::LocallyActiveDataCoarsenSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > crse_level_mgr,
      tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > fine_level_mgr,
      xfer::LocallyActiveDataCoarsenPatchStrategy<DIM>* patch_strategy =
            ((xfer::LocallyActiveDataCoarsenPatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::LocallyActiveDataCoarsenTransactionFactory<DIM> > 
            transaction_factory = 
            (xfer::LocallyActiveDataCoarsenTransactionFactory<DIM>*)NULL);

   /*!
    * Print the coarsen algorithm state to the specified data stream.
    *
    * @param stream Output data stream. 
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   // the following two functions are not implemented:
   LocallyActiveDataCoarsenAlgorithm(const LocallyActiveDataCoarsenAlgorithm<DIM>&); 
   void operator=(const LocallyActiveDataCoarsenAlgorithm<DIM>&);

   tbox::Pointer< xfer::CoarsenClasses<DIM> > d_coarsen_classes;
   bool d_fill_coarse_data;
   bool d_schedule_created;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataCoarsenAlgorithm.C"
#endif
