//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/RefineAlgorithm.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Refine algorithm for data transfer between AMR levels
//
 
#ifndef included_xfer_RefineAlgorithm
#define included_xfer_RefineAlgorithm

#include "SAMRAI_config.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "RefineClasses.h"
#include "RefineOperator.h"
#include "TimeInterpolateOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RefineTransactionFactory.h"

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/*!
 * @brief Class RefineAlgorithm<DIM> encapsulates the AMR communication
 * pattern to refine data to, copy data to, or fill physical boundary data on 
 * any destination patch level.  The basic procedure for moving data follows 
 * three steps: 
 * -# interpolate data (spatial and possibly temporal) from coarser levels
 * -# copy data from the same level of refinement
 * -# fill physical boundary conditions regions
 * 
 * Each data communication procedure generally consists of three parts: an 
 * algorithm, a schedule, and a patch strategy.  The algorithm describes 
 * the patch data components and time and space interpolation operations,
 * but is independent of the configuration of the patches in an AMR hierarchy.  
 * Patch data items and their associated spatial and time interpolation 
 * operators are registered with an instantiation of this algorithm class.  
 * To generate the communication dependencies for a particular patch hierarchy 
 * configuration, the algorithm creates a refine schedule based on the state 
 * of a given hierarchy and the information in the algorithm.  The schedule 
 * can then perform the communication operations that move data to the
 * destination patch level using the associated operators.  User-defined
 * operations (such as filling physical boundaries and special interpolation
 * procedures) are provided through a refine patch strategy object.
 *
 * In general, source data is copied into the designated scratch data for
 * temporary processing.  The scratch space must contain sufficient ghost
 * cells to accommodate the stencil width of the given interpolation operators
 * and any physical boundary data that must be filled.  The scratch storage is 
 * copied into the destination data space at the end of the communication process.  
 * Thus, copy operations between source, scratch, and destination patch data 
 * objects must be defined.  In general, it is the user's responsibility to 
 * register valid operations with the refine algorithm so that the data 
 * communication can occur.
 *
 * In general, the destination and scratch data components may be the same 
 * (assuming that the scratch component has a sufficient ghost cells width).
 * The source and scratch components SHOULD NOT be the same, since the
 * interiors of the source space would be changed by the use of the scratch
 * data as temporary work space.
 *
 * Note that each refine schedule created by a refine algorithm remains valid
 * as long as the patches involved in the communication process do not change;
 * thus, they can be used for multiple data communication cycles.
 *
 * Typical usage of a refine algorithm to perform inter-patch communication
 * on an AMR hierarchy involves four steps:
 *
 * -# Construct a refine algorithm object.
 * -# Register refine operations with the refine algorithm.  Using the
 *       registerRefine() methods(s), one provides source and destination
 *       patch data information, as well as time and space interpolation
 *       operators as needed.  Two registerRefine() methods appear in this
 *       class; one supports time interpolation, one does not.
 * -# After all operations are registered with the algorithm, one
 *       creates a communication schedule using one of the createSchedule()
 *       methods.  These methods are distinguished by the resulting data
 *       communication pattern (e.g., interpatch communication on a single
 *       level, between two different levels at the same grid resolution, 
 *       interpolation of data between different AMR hierarchy levels, etc.)
 *       Note that when creating a communication schedule, a concrete
 *       instance of a RefinePatchStrategy<DIM> object may be required
 *       to supply physical boundary conditions as well as user-defined
 *       spatial data interpolation operations.
 * -# Invoke the fillData() method in the communication schedule to 
 *       perform the data transfers.        
 *
 * User-defined interpolation operations can be written using the interfaces
 * in RefinePatchStrategy for preprocessRefine() and postProcessRefine().
 * Users who use create these operations must note that all data that is to
 * be used in such operations should be registered with the RefineAlgorithm
 * using registerRefine(), whether or not the data is to be refined.
 *
 * @see xfer::RefineSchedule
 * @see xfer::RefinePatchStrategy
 * @see xfer::RefineClasses
 */
 
template<int DIM> class RefineAlgorithm : public tbox::DescribedClass
{
public:
   /*!
    * Construct a refinement algorithm and initialize its basic state.
    * Note that refinement operations must be registered with this algorithm 
    * before it can do anything useful.  See the registerRefine() routines 
    * for details.
    */
   RefineAlgorithm();
 
   /*!
    * The virtual destructor for the algorithm releases all internal storage.
    */
   virtual ~RefineAlgorithm<DIM>();

   /*!
    * Register a refine operation with the refine algorithm object.  This 
    * routine does not support time interpolation.  Data values will be moved 
    * from the source component to the destination component using scratch
    * component as a temporary work space.  The scratch component must have 
    * sufficient ghost cells to cover the required operator stencil width and 
    * any needed physical boundary ghost cells.  
    * 
    * @param dst       Integer destination patch data index to be filled on the 
    *                  destination level.
    * @param src       Integer source patch data index on the source level.
    * @param scratch   Integer patch data index used as a temporary work space.
    * @param oprefine  tbox::Pointer to refinement operator.  This may be a
    *                  null pointer.  In this case, refinement must be
    *                  handled by the refine patch strategy member functions.
    *                  See the comments for 
    *                  RefinePatchStrategy<DIM>::preprocessRefine() and
    *                  RefinePatchStrategy<DIM>::postprocessRefine().
    */
   void registerRefine(
      const int dst,
      const int src,
      const int scratch,
      tbox::Pointer< RefineOperator<DIM> > oprefine);

   /*!
    * Register a refine operation with the refine algorithm object.  This
    * routine supports time interpolation.  Time interpolation will take place
    * between the old and new source data components on coarser levels.  On the
    * destination level, data will be moved from the source component to the
    * destination component using scratch component as a temporary work space.
    * The scratch component must have sufficient ghost cells to cover the
    * required operator stencil width and any needed physical boundary ghost
    * cells.  The time interpolation operator cannot be null.  When assertion
    * checking is active, passing in a null pointer will result in an
    * unrecoverable exception.
    *
    * @param dst       Integer destination patch data index to be filled on the
    *                  destination level.
    * @param src       Integer source patch data index on the source level.
    * @param src_told  Integer source patch data index for old data used in time interpolation.
    * @param src_tnew  Integer source patch data index for new data used in time interpolation.
    * @param scratch   Integer patch data index used as a temporary work space.
    * @param oprefine  tbox::Pointer to refinement operator.  This may be a
    *                  null pointer.  In this case, refinement must be
    *                  handled by the refine patch strategy member functions.
    *                  See the comments for
    *                  RefinePatchStrategy<DIM>::preprocessRefine() and
    *                  RefinePatchStrategy<DIM>::postprocessRefine().
    * @param optime    tbox::Pointer to time interpolation operator.  This
    *                  pointer may not be null.
    */
   void registerRefine(
      const int dst,
      const int src,
      const int src_told,
      const int src_tnew,
      const int scratch,
      tbox::Pointer< RefineOperator<DIM> > oprefine,
      tbox::Pointer< TimeInterpolateOperator<DIM> > optime);

   /*!
    * @brief Create a communication schedule that moves data from the interiors
    * of the source data components into the interior and boundary cells of the 
    * destination data components on the same level where those sources and
    * destinations overlap.  
    * 
    * Neither time nor spatial interpolation is performed.  
    *
    * Note that the schedule remains valid as long as the patches on the level 
    * do not change; thus, it can be used for multiple data communication
    * cycles.
    *
    * @return tbox::Pointer to refine schedule that performs the data transfers.
    *
    * @param level          tbox::Pointer to level on which interpatch
    *                       transfers occur.  This pointer cannot be null.
    * @param patch_strategy Optional tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       is performed.
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state), 
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
         ((xfer::RefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory = 
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
         ((xfer::RefinePatchStrategy<DIM>*)NULL),
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory = 
         (xfer::RefineTransactionFactory<DIM>*)NULL);


   /*!
    * @brief Create a communication schedule that moves data from the interiors
    * of the source data components on a source level into the interior and 
    * boundary cells of the destination data components on a destination level
    * where those sources and destinations overlap.  
    *
    * Note that both levels must reside in the same AMR hierarchy index space,
    * or in index spaces that represent the same level of mesh refinement.
    * No spatial interpolation is performed.
    *
    * In certain rare cases it may be desired to use this schedule to perform
    * time interpolation, in which case the use_time_interpolation optional
    * argument should be set to true.
    *
    * Note that the schedule remains valid as long as the patches on the levels 
    * do not change; thus, it can be used for multiple data communication
    * cycles.
    * 
    * @return tbox::Pointer to refine schedule that performs the data transfers.
    *
    * @param dst_level      tbox::Pointer to destination level; cannot be null.
    * @param src_level      tbox::Pointer to source level; cannot be null.
    * @param patch_strategy tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       is performed.
    * @param use_time_interpolation Optional boolean flag to create the schedule with
    *                            the ability to perform time interpolation. 
    *                            Default is no time interpolation (false).
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*
    * @brief Similar to the above, except with fill_pattern specified.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Create a communication schedule that moves data from the interiors
    * of the source data components on the patch level and coarser levels in the
    * patch hierarchy into the interior and boundary cells of the destination data
    * components on the given patch level where those sources and destinations overlap.  
    * Data is time interpolated between old and new sources on coarser levels
    * when and where time interpolation is needed and copied from the source
    * components on the patch level into the destination components otherwise.
    *
    * In certain rare cases in may be necessary to perform time interpolation
    * between old and new sources on the given patch level.  In this case
    * the optional argument use_time_interpolation should be set to true.
    * Regardless of the value of this argument, time interpolation on
    * coarser levels will always occur whenever needed.
    *
    * Note that the next coarser level number must correspond to a level
    * in the hierarchy that represents a region of coarser index space than
    * the destination level.
    *
    * Note that the schedule remains valid as long as the patches on the levels 
    * involved in its creation do not change; thus, it can be used for multiple 
    * data communication cycles.
    *
    * @return tbox::Pointer to refine schedule that performs the data transfers.
    *
    * @param level          tbox::Pointer to destination level; cannot be null.
    * @param next_coarser_level Integer number of next coarser patch level in
    *                           the patch hierarchy relative to the destination 
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy      tbox::Pointer to patch hierarchy from which data to
    *                       fill level should come.  This pointer may be null
    *                       only when the next_coarser_level is < 0.  
    * @param patch_strategy tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations and user-defined spatial interpolation
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       or user-defined interpolation is performed.  Note
    *                       that this may cause problems since interpolation of
    *                       data from coarser levels to some finer level may
    *                       require physical boundary data.
    * @param use_time_interpolation Optional boolean flag to create the schedule with
    *                               the ability to perform time interpolation
    *                               Default is no time interpolation (false).
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Create a communication schedule that moves data from the interiors
    * of the source data components on the source level and coarser levels in the 
    * hierarchy into the interior and boundary cells of the destination data 
    * components on the destination level where those sources and destinations 
    * overlap.  Data is time interpolated between old and new sources on coarser 
    * levels when and where time interpolation is needed and from the source data 
    * components on the source level into the destination data components otherwise.
    *
    * This form of schedule construction is typically used after regridding
    * (where the source level is the patch level being replaced by the 
    * destination level in the patch hierarchy) or when the data on destination 
    * patch level is to be overwritten by data interpolated from coarser levels 
    * in the patch hierarchy.  In the first case, data on the destination level
    * will be copied from the source level in regions where those two levels
    * overlap and filled with interpolated values from the hierarchy elsewhere.
    * In the latter case, the source level pointer may be null.  Then, data on 
    * the destination level will be filled using interpolated data from coarser 
    * hierarchy levels.
    *
    * In certain cases it may be necessary to perform time interpolation
    * between old and new sources onto the destination level.  In this case
    * the optional argument use_time_interpolation should be set to true.
    * Regardless of the value of this argument, time interpolation on
    * coarser levels will always occur whenever needed.
    *
    * Note that when the source level pointer is non-null, the index spaces of
    * the source and destination levels must be aligned with one another.
    *
    * Note that the schedule remains valid as long as the patches on the levels 
    * involved in its creation do not change; thus, it can be used for multiple 
    * data communication cycles.
    *
    * @return tbox::Pointer to refine schedule that performs the data transfers.
    *
    * @param dst_level          tbox::Pointer to destination level; cannot be
    *                           null.
    * @param src_level          tbox::Pointer to source level. This pointer may
    *                           be null.  In this case, data on the destination
    *                           level will be filled only using interpolated
    *                           data from coarser hierarchy levels.  When this
    *                           pointer is not null, the source level must live
    *                           in the same AMR hierarchy index space as the
    *                           destination level.
    * @param next_coarser_level Integer number of next coarser patch level in
    *                           a patch hierarchy relative to the destination
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy      tbox::Pointer to patch hierarchy from which data to
    *                       fill level should come.  This pointer may be null
    *                       only when the next_coarser_level is < 0.
    * @param patch_strategy tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations and user-defined spatial interpolation
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       or user-defined interpolation is performed.  Note
    *                       that this may cause problems since interpolation of
    *                       data from coarser levels to some finer level may
    *                       require physical boundary data.
    * @param use_time_interpolation Optional boolean flag to create the schedule with
    *                               the ability to perform time interpolation
    *                               Default is no time interpolation (false).
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Similar to the version of createSchedule without the
    * @c fill_pattern argument.
    * 
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    * 
    * @param level          tbox::Pointer to destination level; cannot be
    *                           null.
    * @param next_coarser_level Integer number of next coarser patch level in
    *                           a patch hierarchy relative to the destination
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy      tbox::Pointer to patch hierarchy from which data to
    *                       fill level should come.  This pointer may be null
    *                       only when the next_coarser_level is < 0.
    * @param patch_strategy tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations and user-defined spatial interpolation
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       or user-defined interpolation is performed.  Note
    *                       that this may cause problems since interpolation of
    *                       data from coarser levels to some finer level may
    *                       require physical boundary data.
    * @param use_time_interpolation Optional boolean flag to create the schedule with
    *                               the ability to perform time interpolation
    *                               Default is no time interpolation (false).
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    *
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Similar to the version of createSchedule without the
    * @c fill_pattern argument.
    *
    * @param fill_pattern Indicates which parts of the destination level
    * to fill.  See RefineSchedule for valid values.
    *
    * @param dst_level          tbox::Pointer to destination level; cannot be
    *                           null.
    * @param src_level          tbox::Pointer to source level. This pointer may
    *                           be null.  In this case, data on the destination
    *                           level will be filled only using interpolated
    *                           data from coarser hierarchy levels.  When this
    *                           pointer is not null, the source level must live
    *                           in the same AMR hierarchy index space as the
    *                           destination level.
    * @param next_coarser_level Integer number of next coarser patch level in
    *                           a patch hierarchy relative to the destination
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy      tbox::Pointer to patch hierarchy from which data to
    *                       fill level should come.  This pointer may be null
    *                       only when the next_coarser_level is < 0.
    * @param patch_strategy tbox::Pointer to a refine patch strategy that
    *                       provides user-defined physical boundary filling
    *                       operations and user-defined spatial interpolation
    *                       operations.  If this patch strategy is null
    *                       (default state), then no physical boundary filling
    *                       or user-defined interpolation is performed.  Note
    *                       that this may cause problems since interpolation of
    *                       data from coarser levels to some finer level may
    *                       require physical boundary data.
    * @param use_time_interpolation Optional boolean flag to create the schedule with
    *                               the ability to perform time interpolation
    *                               Default is no time interpolation (false).
    * @param transaction_factory Optional tbox::Pointer to a refine transaction
    *                            factory that creates data transactions for the
    *                            schedule.  If this pointer is null (default state),
    *                            then a StandardRefineTransactionFactory object
    *                            will be used.
    *
    */
   tbox::Pointer< xfer::RefineSchedule<DIM> > createSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      xfer::RefinePatchStrategy<DIM>* patch_strategy =
            ((xfer::RefinePatchStrategy<DIM>*)NULL),
      bool use_time_interpolation = false,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory =
         (xfer::RefineTransactionFactory<DIM>*)NULL);

   /*!
    * @brief Given a previously-generated refine schedule, check for consistency
    * with this refine algorithm object to see whether a call to resetSchedule()
    * is a valid operation.
    * 
    * @return Boolean true if schedule reset is valid; false otherwise.
    *
    * @param schedule  tbox::Pointer to refine schedule, which cannot be null.
    */
   bool checkConsistency(tbox::Pointer< xfer::RefineSchedule<DIM> > schedule) const;

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
    * data centering, ghost cell widths, stencil requirements, etc.).  However,
    * However, the specific source, destination, scratch patch data ids and
    * operators can be different.  Detailed and fairly complete error checking
    * is performed when this routine is called to prevent potential errors or
    * unexpected behavior. 
    *
    * @param schedule  tbox::Pointer to refine schedule, which cannot be null.
    */
   void resetSchedule(tbox::Pointer< xfer::RefineSchedule<DIM> > schedule) const;

   /*!
    * @brief Return const reference to the pointer to refine equivalence classes
    * used in algorithm.
    */
   const tbox::Pointer< RefineClasses<DIM> >& getEquivalenceClasses() const;

   /*!
    * @brief Set the pointer to the refine equivalence classes to be equal
    * to the given argument.
    *
    * @param refine_classes A pointer to refine equivalence classes
    */
   void setEquivalenceClasses(
      const tbox::Pointer< RefineClasses<DIM> > refine_classes);

   /*!
    * Print the refine algorithm state to the specified data stream.
    *
    * @param stream Output data stream.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   RefineAlgorithm(const RefineAlgorithm<DIM>&);	// not implemented
   void operator=(const RefineAlgorithm<DIM>&);	// not implemented

   tbox::Pointer< RefineClasses<DIM> > d_refine_classes;

   bool d_schedule_created;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefineAlgorithm.C"
#endif
