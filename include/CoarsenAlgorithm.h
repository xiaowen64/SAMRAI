//
// File:	CoarsenAlgorithm.h
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Coarsening algorithm for data transfer between AMR levels
//
 
#ifndef included_xfer_CoarsenAlgorithm
#define included_xfer_CoarsenAlgorithm

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_xfer_CoarsenClasses
#include "CoarsenClasses.h"
#endif
#ifndef included_xfer_CoarsenOperator
#include "CoarsenOperator.h"
#endif
#ifndef included_xfer_CoarsenSchedule
#include "CoarsenSchedule.h"
#endif
#ifndef included_xfer_CoarsenPatchStrategy
#include "CoarsenPatchStrategy.h"
#endif
#ifndef included_xfer_CoarsenTransactionFactory
#include "CoarsenTransactionFactory.h"
#endif

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/*!
 * @brief Class CoarsenAlgorithm<DIM> encapsulates the AMR communication
 * pattern to coarsen data from a finer level to any coarser level.  
 * Most often, data is coarsened from the interiors of source patch 
 * components on the source patch level into interiors of destination 
 * patch components on the destination level.  See comments for the
 * coarsen algorithm constructor for variations that are possible 
 * for (adventurous?) users. If the coarsening operators require ghost cells 
 * on a source component, then sufficient ghost cell storage must be provided 
 * by the source patch data component, and those ghost cells must be filled 
 * before calling the data coarsening routines.
 *
 * Communication algorithms generally consist of three parts: an algorithm,
 * a schedule, and a patch strategy.  The algorithm describes the communication
 * between patch data items but is independent of the configuration of the
 * AMR hierarchy.  hier::Patch data items and their associated coarsening
 * operators are registered with the algorithm.  To generate the communication
 * dependencies for a particular hierarchy configuration, the algorithm
 * generates a schedule based on the current hierarchy configuration.  This 
 * schedule then performs the communication based on the registered data types 
 * and their associated operators.  User-defined pre-processing and
 * post-processing are provided through the abstract patch strategy class.
 *
 * The source patch data space is used during processing to store temporary
 * data.  Thus, the user-defined coarsening operators should operate on the
 * source space by using the patch data with those indices.
 *
 * Note that each coarsen schedule created by a coarsen algorithm remains valid
 * as long as the patches involved in the communication process do not change;
 * thus, they can be used for multiple data communication cycles.
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
 *       CoarsenPatchStrategy<DIM> object may be required to supply 
 *       user-defined spatial data coarsening operations.
 * -# Invoke the coarsenData() method in the communication schedule to
 *       perform the data transfers.
 *
 * User-defined coarsening operations can be written using the interfaces
 * in CoarsenPatchStrategy for preprocessCoarsen() and postProcessCoarsen().
 * Users who use create these operations must note that all data that is to
 * be used in such operations should be registered with the CoarsenAlgorithm
 * using registerCoarsen(), whether or not the data is to be coarsened.
 *
 * @see xfer::CoarsenSchedule
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::CoarsenClasses
 */
 
template<int DIM> class CoarsenAlgorithm  : public tbox::DescribedClass
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
    */
   CoarsenAlgorithm(bool fill_coarse_data = false);
 
   /*!
    * The virtual destructor for the algorithm releases all internal storage.
    */
   virtual ~CoarsenAlgorithm<DIM>();

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
    * between coarse and fine levels is four in this case, then the source patch
    * data is required to have at least eight ghost cells.
    * 
    * @param dst       hier::Patch data index filled on destination level.
    * @param src       hier::Patch data index coarsened from the source level.
    * @param opcoarsen tbox::Pointer to coarsening operator.  This may be a
    *                  null pointer.
    *                  In this case, coarsening must be handled by the coarsen
    *                  patch strategy member functions.  See the comments for 
    *                  CoarsenPatchStrategy<DIM>::preprocessCoarsen() and
    *                  CoarsenPatchStrategy<DIM>::postprocessCoarsen(). 
    * @param gcw_to_coarsen Integer vector ghost cell width when data should 
    *                       be coarsened from ghost cell regions of the source 
    *                       (finer) level into the destination (coarser) level.
    *                       By default, it is the vector of zeros indicating
    *                       that data should be coarsened from from patch
    *                       interiors on the source level.  If this argument is
    *                       used, its value should be the cell width, in the
    *                       destination (coarser) level index space, of the
    *                       region around the fine level where this coarsening
    *                       should occur.  This argument should only be provided
    *                       by those who specifically require this special
    *                       behavior and know how to properly process the patch
    *                       data on coarse and fine patch levels during the 
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
    * do not change; thus, it can be used for multiple data communication
    * cycles.
    *
    * @return tbox::Pointer to coarsen schedule that performs the data
    *         transfers. 
    * 
    * @param crse_level     tbox::Pointer to coarse (destination) level.
    * @param fine_level     tbox::Pointer to fine (source) level.
    * @param patch_strategy tbox::Pointer to a coarsen patch strategy that
    *                       provides user-defined coarsen operations.  If this
    *                       patch strategy is null (default state), then no 
    *                       user-defined coarsen operations will be performed.  
    * @param transaction_factory Optional tbox::Pointer to a coarsen transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardCoarsenTransactionFactory object
    *                            will be used.
    */
   tbox::Pointer< xfer::CoarsenSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
      tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
      xfer::CoarsenPatchStrategy<DIM>* patch_strategy =
            ((xfer::CoarsenPatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::CoarsenTransactionFactory<DIM> > transaction_factory =
            (xfer::CoarsenTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Given a previously-generated coarsen schedule, check for
    * consistency with this coarsen algorithm object to see whether a call to
    * resetSchedule is a valid operation.
    *
    * @return Boolean true if schedule reset is valid; false otherwise.
    *
    * @param schedule  tbox::Pointer to coarsen schedule, which cannot be null.
    */
   bool checkConsistency(
      tbox::Pointer< xfer::CoarsenSchedule<DIM> > schedule) const;

   /*!
    * @brief Given a previously-generated coarsen schedule, reconfigure it to
    * peform the communication operations registered with this coarsen algorithm
    * object.  That is, the schedule will be transformed so that it will
    * function as though this coarsen algorithm created it.  Note that the set
    * of operations registered with this coarsen algorithm must be essentially
    * the same as those registered with the coarsen algorithm that created the
    * schedule originally.
    * That is, the number of operations registered must be the same and the
    * source and destination patch data items and operators must have identical
    * characteristics (i.e., data centering, ghost cell widths, stencil
    * requirements, etc.).  However, the specific source, destination patch data
    * ids and operators can be different.  Detailed and fairly complete error
    * checking is performed when this routine is called to prevent potential
    * errors or unexpected behavior.
    *
    * @param schedule  tbox::Pointer to coaren schedule, which cannot be null.
    */
   void resetSchedule(
      tbox::Pointer< xfer::CoarsenSchedule<DIM> > schedule) const;

   /*!
    * Print the coarsen algorithm state to the specified data stream.
    *
    * @param stream Output data stream. 
    */
   void printClassData(ostream& stream) const;

private:
   CoarsenAlgorithm(const CoarsenAlgorithm<DIM>&); // not implemented
   void operator=(const CoarsenAlgorithm<DIM>&);	  // not implemented

   tbox::Pointer< xfer::CoarsenClasses<DIM> > d_coarsen_classes;
   bool d_fill_coarse_data;
   bool d_schedule_created;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenAlgorithm.C"
#endif
