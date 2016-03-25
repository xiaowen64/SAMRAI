//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataRefineSchedule.h $
// Package:     SAMRAI data transfer
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1818 $
// Modified:	$LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
// Description:	Refine schedule for locally-active data transfer between AMR levels
//
 
#ifndef included_xfer_LocallyActiveDataRefineSchedule
#define included_xfer_LocallyActiveDataRefineSchedule

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_LocallyActiveDataPatchLevelManager
#include "LocallyActiveDataPatchLevelManager.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_GridGeometry
#include "GridGeometry.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Schedule
#include "tbox/Schedule.h"
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif
#ifndef included_xfer_LocallyActiveDataFillBoxSet
#include "LocallyActiveDataFillBoxSet.h"
#endif
#ifndef included_xfer_RefineClasses
#include "RefineClasses.h"
#endif
#ifndef included_xfer_LocallyActiveDataRefinePatchStrategy
#include "LocallyActiveDataRefinePatchStrategy.h"
#endif
#ifndef included_xfer_LocallyActiveDataRefineTransactionFactory
#include "LocallyActiveDataRefineTransactionFactory.h"
#endif

namespace SAMRAI {
   namespace xfer {

/*!
 * @brief Class LocallyActiveDataRefineSchedule<DIM> performs the communication 
 * operations that refine data to, copy data to, or fill physical boundary 
 * data on a destination patch level, where the data may be defined only on some
 * patches (i.e., the patch data is "locally-active").  This class is based on the 
 * RefineSchedule<DIM> class.  However, it has a reduced set of functionality since
 * it treats locally-active data.  For example, this class does not support time
 * interpolation, and there is only one version of the constructor (vs. three for 
 * RefineSchedule<DIM>), the schedule cannot be reset, etc.
 *
 * Source data is copied into the provided scratch space for temporary
 * processing.  The scratch space must contain sufficient ghost
 * cells to accommodate the stencil width of the given interpolation operators
 * and any physical boundary data that must be filled.  The scratch data
 * is copied into the destination space at the end of the process.
 * The communication schedule is executed by calling member function fillData().
 *
 * Each schedule object is typically created by a refine algorithm and
 * represents the communication dependencies for a particular configuration
 * of the AMR hierarchy.  The communication schedule is only valid for that
 * particular configuration and must be regenerated when the AMR patch
 * hierarchy changes.  However, as long as the patch levels involved in
 * the creation of the schedule remain unchanged, the schedule may be used
 * for multiple communication cycles.  For more information about creating
 * refine schedules for locally-active patch data, see the 
 * LocallyActiveDataRefineAlgorithm<DIM> header file.
 *
 * NOTE: Algorithmic variations are available by calling the static method
 *       LocallyActiveDataRefineSchedule<DIM>::setScheduleGenerationMethod(),
 *       which sets the option for all instances of the class.
 *
 * @see xfer::LocallyActiveDataRefineAlgorithm
 * @see xfer::LocallyActiveRefinePatchStrategy
 * @see xfer::RefineClasses
 */

template<int DIM> class LocallyActiveDataRefineSchedule :
public tbox::DescribedClass
{
public:

   /*!
    * Static function to set box intersection algorithm to use during
    * schedule construction for all LocallyActiveDataRefineSchedule<DIM> 
    * objects. If this method is not called, the default will be used.
    *
    * Note that the ability to change the method from the default case is
    * disabled currently.
    *
    * @param  method   string identifying box intersection method.  Valid
    *                  choices are:  "BOX_TREE" (default case), "BOX_GRAPH",
    *                  and "ORIG_NSQUARED".   More details can be found below
    *                  in the comments for the generateSchedule() routine.
    *
    * If an invalid string is passed, an unrecoverable error will result.
    */
   static void setScheduleGenerationMethod(const std::string& method);

   /*!
    * @brief Constructor to create a refine schedule that moves data from the
    * interiors of source patch data components on the source level into
    * the interiors and ghost cells of destination patch data components
    * on the destination level.
    *
    * Only data on the intersection of the source and destination patch components
    * will be copied.  The source and destination patch levels must reside in the
    * same index space.  However, the levels do not have to be in the same AMR patch
    * hierarchy.  Generally, this constructor is called by a 
    * LocallyActiveDataRefineAlgorithm<DIM> object.
    *
    * @param dst_level       Pointer to destination patch level.
    * @param dst_level_mgr   Pointer to destination level data manager; cannot be null.
    * @param src_level       Pointer to source patch level.
    * @param src_level_mgr   Pointer to source level data manager; cannot be null.
    * @param refine_classes  Pointer to structure containing patch data and
    *                        operator information.  In general, this is
    *                        constructed by the calling RefineAlgorithm<DIM>
    *                        object.  This pointer cannot be null.
    * @param transaction_factory  Pointer to a factory object that will create
    *                        data transactions; cannot be null.
    * @param patch_strategy  Pointer to a refine patch strategy object that
    *                        provides user-defined physical boundary filling
    *                        operations.   This pointer may be null, in which
    *                        case no boundary filling operations will occur.
    * @param use_time_interpolation Optional boolean flag indicating whether to use time
    *                        interpolation when setting data on the destination level.
    *                        Default is no time interpolation.
    *
    * When assertion checking is active, an unrecoverable assertion will result
    * if a null pointer is passed as described above.
    */
   LocallyActiveDataRefineSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory,
      xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
      bool use_time_interpolation = false);

   /*!
    * @brief Constructor to create a refine schedule that moves data from the
    * interiors of source patch data components on the source level and
    * coarser levels in the patch hierarchy into the interiors and ghost
    * cells of destination patch data components on the destination level.
    *
    * Only data on the intersection of the source and destination patch
    * components, where both data are defined, will be copied.  If portions
    * of the destination level remain unfilled, then the algorithm recursively
    * fills those unfilled portions from coarser levels in the AMR hierarchy.
    * The source and destination patch levels must reside in the same index space.
    * However, the levels do not have to be in the same AMR patch hierarchy.
    * In general, this constructor is called by a
    * LocallyActiveDataRefineAlgorithm<DIM> object.
    *
    * @param dst_level       Pointer to destination patch level; cannot be null.
    * @param dst_level_mgr   Pointer to destination level data manager; cannot be null.
    * @param src_level       Pointer to source patch level; must be in same
    *                        index space as destination level.  This pointer
    *                        may be null, in which case the destination level
    *                        will be filled only using data interpolated from
    *                        coarser levels in the AMR hierarchy.
    * @param src_level_mgr   Pointer to source level data manager; may be null
    *                        only if src_level pointer is null.
    * @param next_coarser_level Integer number of next coarser level in
    *                           AMR patch hierarchy relative to the destination
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy       Pointer to patch hierarchy.  This pointer may be
    *                        null only if the next_coarser_level number is < 0,
    *                        indicating that there is no level in the hierarchy
    *                        coarser than the destination level.
    * @param refine_classes  Pointer to structure containing patch data and
    *                        operator information.  In general, this is constructed 
    *                        by the calling LocallyActiveDataRefineAlgorithm<DIM>
    *                        object.  This pointer cannot be null.
    * @param transaction_factory  Pointer to a factory object that will create
    *                        data transactions; cannot be null.
    * @param patch_strategy  Pointer to a refine patch strategy object that
    *                        provides user-defined physical boundary filling
    *                        operations.   This pinter may be null, in which
    *                        case no boundary filling or user-defined refine
    *                        operations will occur.
    * @param use_time_interpolation Optional boolean flag indicating whether to use time
    *                        interpolation when setting data on the destination level.
    *                        Default is no time interpolation.
    *
    * When assertion checking is active, an unrecoverable assertion will result
    * if a null pointer is passed as described above.
    */
   LocallyActiveDataRefineSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory,
      xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
      bool use_time_interpolation = false);

   /*!
    * Virtual destructor for the schedule releases all internal storage.
    */
   virtual ~LocallyActiveDataRefineSchedule<DIM>();

   /*!
    * @brief Execute the stored communication schedule and perform
    * the data movement.
    *
    * @param fill_time                 Double time for filling operation.
    * @param do_physical_boundary_fill Boolean flag used to bypass the
    *                                  physical boundary data filling
    *                                  operations on the destination level.
    *                                  The default value is true indicating
    *                                  that boundary data will be filled
    *                                  (assuming a non-null refine patch
    *                                  strategy pointer was passed to the
    *                                  createSchedule() function.  Note that
    *                                  even when the value is false, boundary
    *                                  routines may be called on levels coarser
    *                                  than the destination level if such data
    *                                  is needed for proper interpolation.
    */
   void fillData(double fill_time,
                 bool do_physical_boundary_fill = true) const;

   /*!
    * Print the refine schedule data to the specified data stream.
    *
    * @param stream Output data stream.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   /*
    * The following two functions are not implemented. 
    */
   LocallyActiveDataRefineSchedule(const LocallyActiveDataRefineSchedule&);
   void operator=(const LocallyActiveDataRefineSchedule&);

   /*!
    * @brief This private constructor creates a communication schedule
    * that fills the destination level on the specified fill boxes only.
    *
    * This constructor is used by the refine schedule algorithm during the
    * recursive schedule generation process.  Internal flags are
    * to reflect that fact.
    */
   LocallyActiveDataRefineSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory,
      const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
      xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy);

   /*!
    * @brief Finish the schedule construction by recursing to coarser levels as needed.
    */
   void finishScheduleConstruction(
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
      bool use_time_interpolation);

   /*!
    * @brief Recursively fill the destination level with data at the
    * given time.
    *
    * @param fill_time solution time when the fill takes place
    * @param do_physical_boundary_fill Boolean indicating whether
    * to call user-supplied boundary filling
    * routines regardless of whether this is needed based on ghost cell
    * width of destination data components or stencil width of some
    * interpolation operator.
    */
   void recursiveFill(double fill_time,
                      bool do_physical_boundary_fill) const;

   /*!
    * @brief Fill the physical boundaries for each patch on the
    * specified patch level.
    *
    * @param level     level to fill physical boundaries
    * @param fill_time solution time when the fill takes place
    */
   void fillPhysicalBoundaries(tbox::Pointer< hier::PatchLevel<DIM> > level,
                               double fill_time) const;

   /*!
    * @brief Allocate scratch space on the specified level and
    * return the allocated patch data indices in the allocate
    * manager for later deallocation.
    */
   void allocateScratchSpace(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
      double fill_time,
      hier::LocallyActiveDataPatchLevelManager<DIM>& allocate_mgr) const;

   /*!
    * @brief Copy the scratch space into the destination space.
    *
    * If the scratch and destination spaces are the same,
    * then no copying is performed.
    */
   void copyScratchToDestination(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr) const;

   /*!
    * @brief Refine scratch data between coarse and fine patch levels.
    */
   void refineScratchData() const;

   /*
    * Generate a communication schedule between the source and destination
    * patch levels.  This algorithm only fills into fill boxes on the
    * destination level.  The boxes that remain unfilled are returned.
    */
   void generateCommunicationSchedule(
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
      const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
      tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_unfilled_boxes,
      const bool use_time_interpolation);

   /*!
    * @brief Calculate the default fill boxes for the specified patch level.
    *
    * The default fill boxes cover the interiors plus the ghost cells.
    */
   void allocateDefaultFillBoxes(
      tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
      const hier::IntVector<DIM>& fill_ghost_width);

   /*!
    * @brief Calculate the maximum ghost cell width of all destination
    * patch data components.
    */
   hier::IntVector<DIM> getMaxDestinationGhosts() const;

   /*!
    * @brief Calculate the maximum ghost cell width of all scratch patch data
    * components.
    */
   hier::IntVector<DIM> getMaxScratchGhosts() const;

   /*!
    * @brief Calculate the maximum ghost cell width required for all stencils.
    */
   hier::IntVector<DIM> getMaxStencilGhosts() const;

   /*!
    * @brief This function is called from each constructor to cache local copies
    * of hierachy information and to compute the necessary scratch data,
    * destination data, and interpolation stencil ghost cell widths used
    * during schedule construction.
    */
   void initializeDomainAndGhostInformation(
      bool recursive_schedule);

   /*!
    * @brief Utility function to set up local copies of patch data source,
    * destination, etc. indices and necessary data interpolation information
    * stored in the refine classes object generated by the refine algorithm.
    *
    * An array of refine data items is stored locally here to facilitate
    * interaction with transations.
    */
   void setRefineItems(const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes);

   /*!
    * @brief Utility function to clear local copies of refine items.
    */
   void clearRefineItems();

   /*!
    * @brief Utility function to check whether scratch data items in collection
    * of refine classes have sufficient ghost cell widths to handle
    * user-defined interpolation operations.
    *
    * Although it is redundant if the refine algorithm created the refine classes,
    * the routine xfer::RefineClasses<DIM>::checkRefineItem() is also called.
    *
    * If any entries are erroneous an assertion is thrown with a descriptive
    * error message and program halts.
    */
   void initialCheckRefineClassItems() const;

   /*!
    * Constant int vectors used to avoid recreating these vectors in loops.
    */
   static const hier::IntVector<DIM> s_constant_zero_intvector;
   static const hier::IntVector<DIM> s_constant_one_intvector;

   /*!
    * Selects algorithm used to generate communication schedule.
    */
   static std::string s_schedule_generation_method;

   static bool s_printing;

   /*!
    * Structures that store refine data items.
    */
   tbox::Pointer< xfer::RefineClasses<DIM> > d_refine_classes;
   int d_number_refine_items;
   const typename xfer::RefineClasses<DIM>::Data** d_refine_items;

   /*!
    * Cached pointer to the destination patch level and data manager.  
    */
   tbox::Pointer< hier::PatchLevel<DIM> > d_dst_level;
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > d_dst_level_mgr;

   /*!
    * Object supporting interface to user-defined boundary filling and
    * spatial data interpolation operations.
    */
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* d_refine_patch_strategy;

   /*!
    * Factory object used to create data transactions when schedule is constructed.
    */
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > d_transaction_factory;

   /*!
     Cached copy of maximum stencil ghost cell widths.
    */
   hier::IntVector<DIM> d_max_stencil_gcw;

   /*!
    * Cached copy of maximum scratch ghost cell widths.
    */
   hier::IntVector<DIM> d_max_scratch_gcw;

   /*!
    * Width of ghost cell region to fill passed to user supplied
    * physical boundary condition routine.
    */
   hier::IntVector<DIM> d_boundary_fill_ghost_width;

   /*!
    * Flag indicating whether user's physical boundary data filling
    * routine should be forced at last step of level filling process.
    *
    * This flag is true when doing recursive filling, because the ghost
    * data may be needed by finer levels (regardless of whether the user
    * requested ghost boundary filling).  This variable is set in
    * the constructors, which knows whether the object is being constructed
    * for recursive filling.
    *
    * For efficiency, we only force boundary filling when, during object
    * construction, we determine that the ghost cells do exist.
    */
   bool d_force_boundary_fill;

   /*!
    * Boolean flag indicating whether physical domain
    * is comprised as a single box region.
    */
   bool d_domain_is_one_box;

   /*!
    * Cached box describing physical domain when that domain
    * is comprised as a single box region
    */
   hier::Box<DIM> d_domain_box;

   /*!
    * Number of non-zero entries in periodic shift vector.
    */
   int d_num_periodic_directions;

   /*!
    * Cached copy of the periodic shift vector.
    */
   hier::IntVector<DIM> d_periodic_shift;

   /*!
    * Level-to-level communication schedule between the source and destination.
    *
    * d_coarse_priority_level_schedule handles
    * the situation where coarse data should take precedence at
    * coarse-fine boundaries for data types holding values at patch
    * boundaries but which are considered interior values.
    * d_fine_priority_level_schedule handles the situation where
    * fine data should take precedence.
    */
   tbox::Pointer<tbox::Schedule> d_coarse_priority_level_schedule;

   /*!
    * Level-to-level communication schedule between the source and destination.
    *
    * d_coarse_priority_level_schedule handles
    * the situation where coarse data should take precedence at
    * coarse-fine boundaries for data types holding values at patch
    * boundaries but which are considered interior values.
    * d_fine_priority_level_schedule handles the situation where
    * fine data should take precedence.
    */
   tbox::Pointer<tbox::Schedule> d_fine_priority_level_schedule;

   /*!
    * Schedule to recursively fill data from the next coarser hierarchy level.
    *
    * This schedule describes how to fill the coarser level so that the coarse
    * data can be interpolated into the fine fill boxes on the destination.
    * If no coarser data is needed to fill the fill boxes on the destination
    * level, then this pointer is NULL.
    */
   tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> > d_coarse_schedule;

   /*!
    * Pointer to temporary coarser level necessary to interpolate
    * data into the fill boxes of the destination.
    *
    * This coarser level is filled by the refine schedule above.  If no coarser
    * level data is needed, then this pointer will be NULL.  Note that the coarser
    * level may not have the same mapping as the destination level; see the mapping
    * array below.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > d_coarse_level;

   /*!
    * Coarse level active data manager for treating locally-active
    * data on temporary coarse level.  If no coarser level data is needed, then 
    * the manager pointer will be NULL.
    */
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > d_coarse_level_mgr;

   /*!
    * The mapping between the patches on the coarser level
    * and the patches on the fine level.
    *
    * For each destination patch on the fine level, there
    * may be zero or more coarse patches, depending on the fill boxes for
    * the destination patch and also any extending to the physical boundaries
    * required to maintain proper box relationships.
    */
   tbox::Array<int> d_coarse_to_fine_mapping;

   /*!
    * An array of fill boxes used when refining between coarse and fine patches.
    *
    * This array will have one entry for each local patch on the coarse patch level.
    */
   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> > d_la_fine_fill_boxes;

   /*!
    * Arrays for overlaps and source mask boxes used in construction of transactions.
    *
    * They are declared in the class to make memory management efficient
    * since that function is called many times.
    *
    * The size of these arrays is controlled by d_max_fill_boxes, which is
    * is set by taking the max over the number of fill boxes for each
    * destination patch to be filled.
    */
   tbox::Array< tbox::Pointer< hier::BoxOverlap<DIM> > > d_overlaps;
   int d_max_fill_boxes;
   hier::BoxArray<DIM> d_src_masks;

   /*!
    * Timer objects for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_fill_data;
   tbox::Pointer<tbox::Timer> t_gen_comm_sched;
   tbox::Pointer<tbox::Timer> t_finish_sched_const;

   tbox::Pointer<tbox::Timer> t_gen_comm_sched_unfilled;
   tbox::Pointer<tbox::Timer> t_gen_comm_sched_trans;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataRefineSchedule.C"
#endif
