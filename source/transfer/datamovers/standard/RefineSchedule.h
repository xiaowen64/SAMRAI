//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/RefineSchedule.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Refine schedule for data transfer between AMR levels
//
 
#ifndef included_xfer_RefineSchedule
#define included_xfer_RefineSchedule

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "BoxList.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Schedule.h"
#include "tbox/Timer.h"
#include "FillBoxSet.h"
#include "RefineClasses.h"
#include "RefinePatchStrategy.h"
#include "RefineTransactionFactory.h"

namespace SAMRAI {
   namespace xfer {

/*!
 * @brief Class RefineSchedule<DIM> performs the communication
 * operations that refine data to, copy data to, or fill physical boundary
 * data on a destination patch level.
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
 * refine schedules, see the RefineAlgorithm<DIM> header file.
 *
 * NOTE: Algorithmic variations are available by calling the static method
 *       RefineSchedule<DIM>::setScheduleGenerationMethod(), which
 *       sets the option for all instances of the class.
 *
 * Some constructors accept the argument @c fill_pattern.  This
 * string controls which types of cells are filled and which are
 * omitted from the filling process.  Valid values are:
 * - @c "DEFAULT_FILL" Fill interior and ghost cells.
 * - @c "FILL_INTERIORS_ONLY" Fill interior cells only.
 * - @c "FILL_LEVEL_BORDERS_ONLY" Fill ghosts on level borders only.
 * - @c "FILL_LEVEL_BORDERS_AND_INTERIORS" Fill interior and
 *      ghosts on level borders.
 *
 * @see xfer::RefineAlgorithm
 * @see xfer::RefinePatchStrategy
 * @see xfer::RefineClasses
 */
 
template<int DIM> class RefineSchedule : public tbox::DescribedClass
{
public:
   /*!
    * Static function to set box intersection algorithm to use during
    * schedule construction for all RefineSchedule objects.
    * If this method is not called, the default will be used.
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
    * hierarchy.  Generally, this constructor is called by a RefineAlgorithm<DIM> object.
    *
    * @param fill_pattern    Indicates which parts of the destination level
    *                        to fill.  See RefineSchedule for valid values.
    * @param dst_level       Pointer to destination patch level.
    * @param src_level       Pointer to source patch level.
    * @param refine_classes  Pointer to structure containing patch data and
    *                        operator information.  In general, this is
    *                        constructed by the calling RefineAlgorithm<DIM>
    *                        object.
    * @param transaction_factory  Pointer to a factory object that will create
    *                        data transactions.
    * @param patch_strategy  Pointer to a refine patch strategy object that
    *                        provides user-defined physical boundary filling
    *                        operations.   This pointer may be null, in which 
    *                        case no boundary filling operations will occur.
    * @param use_time_interpolation Optional boolean flag indicating whether to use time
    *                        interpolation when setting data on the destination level.
    *                        Default is no time interpolation.
    *
    * When assertion checking is active, unrecoverable assertions will result
    * if either patch level pointer, the refine classes pointer, or the transaction
    * factory pointer, is null.
    */
   RefineSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
      xfer::RefinePatchStrategy<DIM>* patch_strategy,
      bool use_time_interpolation = false);

   /*!
    * @brief Constructor to create a refine schedule that moves data from the
    * interiors of source patch data components on the source level and
    * coarser levels in the patch hierarchy into the interiors and ghost 
    * cells of destination patch data components on the destination level.  
    *
    * Only data on the intersection of the source and destination patch 
    * components will be copied.  If portions of the destination level 
    * remain unfilled, then the algorithm recursively fills those unfilled 
    * portions from coarser levels in the AMR hierarchy.  The source
    * and destination patch levels must reside in the same index space.
    * However, the levels do not have to be in the same AMR patch hierarchy.
    * In general, this constructor is called by a RefineAlgorithm<DIM> object.
    *
    * @param fill_pattern    Indicates which parts of the destination level
    *                        to fill.  See RefineSchedule for valid values.
    * @param dst_level       Pointer to destination patch level.
    * @param src_level       Pointer to source patch level; must be in same
    *                        index space as destination level.  This pointer
    *                        may be null, in which case the destination level 
    *                        will be filled only using data interpolated from
    *                        coarser levels in the AMR hierarchy.
    * @param next_coarser_level Integer number of next coarser level in 
    *                           AMR patch hierarchy relative to the destination 
    *                           level.  Note that when the destination level
    *                           has number zero (i.e., the coarsest level), this
    *                           value should be < 0.
    * @param hierarchy       Pointer to patch hierarchy.  This pointer may be
                             null only if the next_coarser_level value is < 0,
    *                        indicating that there is no level in the hierarchy
    *                        coarser than the destination level.
    * @param refine_classes  Pointer to structure containing patch data and
    *                        operator information.  In general, this is
    *                        constructed by the calling RefineAlgorithm<DIM>
    *                        object.
    * @param transaction_factory  Pointer to a factory object that will create
    *                        data transactions.
    * @param patch_strategy  Pointer to a refine patch strategy object that
    *                        provides user-defined physical boundary filling
    *                        operations.   This pinter may be null, in which
    *                        case no boundary filling or user-defined refine
    *                        operations will occur.
    * @param use_time_interpolation Optional boolean flag indicating whether to use time
    *                        interpolation when setting data on the destination level.
    *                        Default is no time interpolation.
    *
    * When assertion checking is active, unrecoverable assertions will result
    * if either patch level pointer, the refine classes pointer, or the transaction
    * factory pointer, is null.
    */
   RefineSchedule(
      const std::string& fill_pattern,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
      xfer::RefinePatchStrategy<DIM>* patch_strategy,
      bool use_time_interpolation = false);

   /*!
    * Virtual destructor for the schedule releases all internal storage.
    */
   virtual ~RefineSchedule<DIM>();

   /*!
    * @brief Reset this refine schedule to perform data transfers
    * asssociated with refine class items in function argument.
    *
    * In general, this function is 
    * called by a RefineAlgorithm<DIM> object.
    *
    * @param refine_classes  Pointer to structure containing patch data and
    *                        operator information.  In general, this is
    *                        constructed by the calling RefineAlgorithm<DIM>
    *                        object.  This pointer must be non-null.
    */
   void reset(const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes);

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
    * @brief Initialize a component selector to store the components
    * needed to allocate source data.
    *
    * @param allocate_vector An empty hier::ComponentSelector that will be
    *                        set to contain the patch data indices for source
    *                        data.
    */
   void initializeSourceVector(hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief Allocate destination space on the destination level and
    * return the allocated patch data indices in the component
    * selector for later deallocation.
    *
    * @param fill_time       Double time for filling operation.
    * @param allocate_vector Component selector that will store the
    *                        allocated patch data indices.
    */
   void allocateDestinationSpace(
      double fill_time,
      hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief Allocate scratch space on the specified level and
    * return the allocated patch data indices in the component
    * selector for later deallocation.
    */
   void allocateScratchSpace(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      double fill_time,
      hier::ComponentSelector& allocate_vector) const;

   /*!
    * Initialize a component selector to store the components needed to
    * allocate source data.
    *
    * @param allocate_vector An empty hier_ComponentSelector that will be
    *                        set to contain the patch data indices for
    *                        destination data.
    */
   void initializeDestinationVector(
      hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief Return const reference to the pointer to
    * refine equivalence classes used in schedule.
    */
   const tbox::Pointer< RefineClasses<DIM> >& getEquivalenceClasses() const;

   /*!
    * Return width of ghost cell region to fill which is passed to user
    * supplied physical boundary condition routine.
    */
   const hier::IntVector<DIM>& getBoundaryFillGhostWidth() const;

   /*!
    * Print the refine schedule data to the specified data stream.
    *
    * @param stream Output data stream.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   RefineSchedule(const RefineSchedule<DIM>&);	// not implemented
   void operator=(const RefineSchedule<DIM>&);		// not implemented

   /*!
    * @brief This private constructor creates a communication schedule
    * that fills the destination level on the specified fill boxes only.
    *
    * This constructor is used by the refine schedule algorithm during the 
    * recursive schedule generation process.  Internal flags are
    * to reflect that fact.
    */
   RefineSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
      tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      xfer::RefinePatchStrategy<DIM>* patch_strategy);

   /*!
    * @brief Finish the schedule construction for the two constructors
    * that take a hierarchy as an argument.
    */
   void finishScheduleConstruction(
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      int next_coarser_level,
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      bool use_time_interpolation,
      bool skip_generate_schedule);

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
    * @param level level to file physical boundaries
    * @param fill_time solution time when the fill takes place
    */
   void fillPhysicalBoundaries(tbox::Pointer< hier::PatchLevel<DIM> > level,
                               double fill_time) const;

   /*!
    * @brief Copy the scratch space into the destination space.
    *
    * If the scratch and destination spaces are the same,
    * then no copying is performed.
    */
   void copyScratchToDestination(tbox::Pointer< hier::PatchLevel<DIM> > level) const;

   /*!
    * @brief Refine scratch data between coarse and fine patch levels.
    */
   void refineScratchData() const;

   /*!
    * @brief Main schedule generation routine which passes control
    * to one of the algorithmic variations based on value of
    * s_default_schedule_generation_method.
    *
    * The resulting transactions will only fill the regions of intersection 
    * between the fill_boxes and the destination level boxes.  The remaining 
    * box regions are returned in unfilled_boxes.
    *
    * The generateCommunicationSchedule() routine invokes various versions 
    * of the schedule generation process implemented in the similarly named 
    * routines below based on the chosen schedule generation method. The 
    * different options will not change the result of the application but 
    * may improve its performance, especially for large numbers of processors.
    * Note that the algorithm choice is set in the getFromInput() method.
    *
    * The possibilities are as follows:
    *
    * - if setScheduleGenerationMethod("BOX_TREE") is called use
    *      generateScheduleBoxTree() to generate the schedule.
    *      NOTE: THIS IS THE DEFAULT OPTION.
    *
    * - if setScheduleGenerationMethod("BOX_GRAPH") is called use
    *      generateScheduleBoxGraph() to generate the schedule.
    *
    * - if setScheduleGenerationMethod("ORIG_NSQUARED") is called use
    *      generateScheduleNSquared() to generate the schedule. 
    */
   void generateCommunicationSchedule(
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
      const bool use_time_interpolation);

   /*!
    * @brief This version of the schedule generation procedure
    * uses N^2 algorithms to determine box intersections;
    * i.e., the original SAMRAI implementation which checks every
    * box against every other.
    */ 
   void generateCommunicationScheduleNSquared(
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
      const bool use_time_interpolation);

   /*!
    * @brief This version of the schedule generation procedure uses a bipartite 
    * graph algorithm to determine which source patches contribute data
    * to each destination patch.
    */
   void generateCommunicationScheduleBoxGraph(
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
      const bool use_time_interpolation);

   /*!
    * @brief This version of the schedule generation procedure uses a recursive
    * binary box tree algorithm to determine which source patches contribute 
    * data to each destination patch and to compute unfilled_boxes.
    */
   void generateCommunicationScheduleBoxTree(
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
      const bool use_time_interpolation);


   /*!
    * @brief Calculate the fill boxes for the specified patch level.
    *
    * Acceptable values for @c fill_pattern is discussed in the
    * class documentation.
    */
   void allocateFillBoxes(
      const std::string& fill_pattern,
      tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Pointer< hier::PatchLevel<DIM> > level,
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
    * @brief Utility function to generate unfilled box regions for a level.
    *
    * This is the original N^2 algorithm in SAMRAI before the development
    * of box tree and box graph algorithms.
    */
   void makeUnfilledBoxesNSquared(
      tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
      const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level) const;

   /*!
    * @brief Function that constructs schedule transactions that 
    * move data from source patch on source level to destination patch
    * on destination level on regions defined by list of fill boxes.
    */
   void constructScheduleTransactions(
      tbox::Pointer<tbox::Schedule> fine_priority_schedule,
      tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
      const hier::BoxList<DIM>& fill_boxes,
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      int dst_patch_id,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      int src_patch_id,
      bool use_time_interpolation);

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

   /*
    * @brief Utility function to clear local copies of refine items.
    */
   void clearRefineItems();

   /*
    * @brief Utility function to check whether scratch data items in collection
    * of refine classes have sufficient ghost cell widths to handle 
    * user-defined interpolation operations.  
    *
    * Although it is redundant if the refine algorithm created the refine classes, 
    * the routine RefineClasses<DIM>::checkRefineItem() is also called.
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

   /*!
    * Structures that store refine data items.
    */
   tbox::Pointer< xfer::RefineClasses<DIM> > d_refine_classes;
   int d_number_refine_items;
   const typename xfer::RefineClasses<DIM>::Data** d_refine_items;

   /*!
    * Cached pointer to the destination patch level.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > d_dst_level;

   /*!
    * Object supporting interface to user-defined boundary filling and
    * spatial data interpolation operations.
    */
   xfer::RefinePatchStrategy<DIM>* d_refine_patch_strategy;

   /*!
    * Factory object used to create data transactions when schedule is constructed.
    */
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > d_transaction_factory;

   /*!
    * Cached copy of maximum stencil ghost cell widths.
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
   tbox::Pointer< xfer::RefineSchedule<DIM> > d_coarse_schedule;

   /*!
    * Pointer to the coarser level necessary to interpolate
    * data into the fill boxes of the destination.
    *
    * This coarser level is filled by the refine schedule above.  If no coarser 
    * level data is needed, then this pointer will be NULL.  Note that the coarser 
    * level may not have the same mapping as the destination level; see the mapping 
    * array below.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > d_coarse_level;

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
    * This array has one entry for each local patch on the coarse patch level.
    */
   tbox::Array< xfer::FillBoxSet<DIM> > d_fine_fill_boxes;

   /*!
    * Arrays for overlaps and source mask boxes used in
    * the private member function constructScheduleTransactions().
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
   tbox::Pointer<tbox::Timer> t_gen_sched_n_squared;
   tbox::Pointer<tbox::Timer> t_gen_sched_box_graph;
   tbox::Pointer<tbox::Timer> t_gen_sched_box_tree; 


};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefineSchedule.C"
#endif
