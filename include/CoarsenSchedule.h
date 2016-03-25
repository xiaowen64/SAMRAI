//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenSchedule.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Coarsening schedule for data transfer between AMR levels
//
 
#ifndef included_xfer_CoarsenSchedule
#define included_xfer_CoarsenSchedule

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "IntVector.h"
#include "ComponentSelector.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Schedule.h"
#include "tbox/Timer.h"
#include "CoarsenClasses.h"
#include "CoarsenPatchStrategy.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "CoarsenTransactionFactory.h"

namespace SAMRAI {
   namespace xfer {

/*!
 * @brief Class CoarsenSchedule<DIM> performs the communication operations
 * to coarsen data from a finer level to a coarser level.  
 * 
 * Typically, data is coarsened from the interiors of source patch components 
 * on the source patch level into interiors of destination patch components on  
 * the destination level.  However, variations are possible for special 
 * situations; see the CoarsenAlgorithm<DIM> class header for more information.
 * Generally, the source patch data must contain sufficient ghost cells to 
 * satisfy the coarsening operators involved.  If a coarsen operator has a 
 * non-zero ghost cell width, then the source ghost cells must be filled before
 * the coarsen schedule is executed.  The communication schedule is executed by 
 * calling member function coarsenData().
 * 
 * Each schedule object is typically created by a coarsen algorithm and
 * represents communication dependencies for a particular configuration
 * of the AMR hierarchy.  The communication schedule is only valid for that
 * particular configuration and must be regenerated when the AMR patch
 * hierarchy changes.  However, as long as the patch levels involved in
 * the creation of the schedule remain unchanged, the schedule may be used
 * for multiple communication cycles.  For more information about creating
 * refine schedules, see the CoarsenAlgorithm<DIM> header file. 
 *
 * NOTE: Algorithmic variations are available by calling the static method
 *       CoarsenSchedule<DIM>::setScheduleGenerationMethod(), which
 *       sets the option for all instances of the class.
 *
 * @see xfer::CoarsenAlgorithm
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::CoarsenClasses
 */
 
template<int DIM> class CoarsenSchedule : public tbox::DescribedClass
{
public:
   /*!
    * Static function to set box intersection algorithm to use during 
    * schedule construction for all CoarsenSchedule objects.  
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
    * Constructor to create a coarsen schedule that coarsens data from 
    * source patch data components on the fine level into the destination patch 
    * data components on the coarse level.  In general, this constructor is 
    * called by a CoarsenAlgorithm<DIM> object.  For possible variations on 
    * data coarsening, see the CoarsenAlgorithm<DIM> class header information.  
    * 
    * If the coarsening operators require data from ghost cells, then the
    * associated source patch data components must have a sufficient ghost
    * cell width and and they must be filled with valid data before calling 
    * coarsenData(). 
    * 
    * @param crse_level        Pointer to coarse (destination) patch level.
    * @param fine_level        Pointer to fine (source) patch level.
    * @param coarsen_classes   Pointer to structure containing patch data and 
    *                          operator information.  In general, this is 
    *                          constructed by the calling CoarsenAlgorithm<DIM>
    *                          object.
    * @param transaction_factory  Pointer to a factory object that will create
    *                          data transactions.
    * @param patch_strategy    Pointer to a coarsen patch strategy object that
    *                          provides user-defined coarsen operations.   This
    *                          ponter may be null, in which case no user-defined 
    *                          coarsen operations will be performed.
    * @param fill_coarse_data  Boolean indicating whether coarse data should
    *                          be filled before coarsening operations are done.
    * 
    * When assertion checking is active, unrecoverable assertions will result
    * if either patch level pointer, the refine classes pointer, or the transaction
    * factory pointer, is null.
    */
   CoarsenSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
      tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
      const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes,
      tbox::Pointer< xfer::CoarsenTransactionFactory<DIM> > transaction_factory,
      xfer::CoarsenPatchStrategy<DIM>* patch_strategy,
      bool fill_coarse_data);

   /*!
    * The virtual destructor for the schedule releases all internal storage.
    */
   virtual ~CoarsenSchedule<DIM>();

   /*!
    * Reset this coarsen schedule to perform data transfers asssociated with
    * coarsen class items in function argument.  In general, this function is
    * called by a CoarsenAlgorithm<DIM> object.
    *
    * @param coarsen_classes  Pointer to structure containing patch data and
    *                         operator information.  In general, this is
    *                         constructed by the calling CoarsenAlgorithm<DIM>
    *                         object.  This pointer must be non-null.
    */
   void reset(const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes);

   /*!
    * Execute the stored communication schedule and perform the data movement.  
    */
   void coarsenData() const;

   /*!
    * Return const reference to the pointer to coarsen equivalence classes
    * used in schedule.
    */
   const tbox::Pointer< xfer::CoarsenClasses<DIM> >& getEquivalenceClasses() const;

   /*!
    * Print the coarsen schedule state to the specified data stream.
    *
    * @param stream Output data stream.
    */
   virtual void printClassData(std::ostream& stream) const;

private:
   CoarsenSchedule(const CoarsenSchedule<DIM>&); // not implemented
   void operator=(const CoarsenSchedule<DIM>&);        // not implemented

   /*!
    * @brief Main schedule generation routine which passes control
    * to one of the algorithmic variations based on value of
    * s_schedule_generation_method. 
    *
    * The resulting communication schedule will move source patch data from a
    * temporary coarse level (i.e., coarsened version of fine level) into the 
    * destination patch data of the destination (coarse) level.
    *
    * The generateSchedule() routine invokes various versions of the schedule 
    * generation process implemented in the similarly named routines below based 
    * on the chosen schedule generation method. The different options will not 
    * change the result of the application but may improve its performance, 
    * especially for large numbers of processors.  Note that the algorithm choice 
    * may be changed by calling the setScheduleGenerationMethod() routine.
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
   void generateSchedule();

   /*!
    * @brief This version of the schedule generation procedure
    * uses N^2 algorithms to determine box intersections;
    * i.e., the original SAMRAI implementation which checks every
    * box against every other.
    */
   void generateScheduleNSquared();

   /*!
    * @brief This version of the schedule generation procedure uses a bipartite
    * graph algorithm to determine which source patches contribute data
    * to each destination patch.
    */
   void generateScheduleBoxGraph();

   /*!
    * @brief This version of the schedule generation procedure uses a recursive
    * binary box tree algorithm to determine which source patches contribute
    * data to each destination patch and to compute unfilled_boxes.
    */
   void generateScheduleBoxTree();

   /*!
    * @brief Generate a temporary coarse level by coarsening the fine level.
    * Note that this function does not allocate patch data storage.
    */
   void generateTemporaryLevel();

   /*!
    * @brief Set up refine algorithms to transfer coarsened data and to fill 
    * temporary coarse level before coarsening operations, if needed.  
    *
    * The associated schedules are set in the generateSchedule() routine.
    */
   void setupRefineAlgorithm();

   /*!
    * @brief Coarsen source patch data from the fine patch level into the 
    * source patch data on the coarse temporary patch level.
    */
   void coarsenSourceData(xfer::CoarsenPatchStrategy<DIM> *patch_strategy) const;

  /*!
    * @brief Calculate the maximum ghost cell width to grow boxes to check
    * for overlaps.
    */
   hier::IntVector<DIM> getMaxGhostsToGrow() const;   

   /*!
    * @brief Function that constructs schedule transactions that
    * move data from source patch on source level to destination patch
    * on destination level.
    */
   void constructScheduleTransactions(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
                                      int dst_patch_id,
                                      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
                                      int src_patch_id);

   /*!
    * @brief Utility function to set up local copies of patch data source, 
    * destination, etc. indices and necessary data coarsening information 
    * stored in the coarsen classes object generated by the coarsen algorithm.  
    * 
    * An array of coarsen data items is stored locally here to facilitate 
    * interaction with transations.
    */
   void setCoarsenItems(const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes);

   /*!
    * @brief Utility function to clear local copies of coarsen items.
    */
   void clearCoarsenItems();

   /*!
    * @brief Utility function to check coarsen items to see whether
    * source and destination patch data components have sufficient ghost
    * cell widths to satisfy the "ghost width to coarsen" functionality
    * described in the CoarsenAlgorithm<DIM> class header.  
    *
    * Specifically, the destination data must have a ghost cell width at least 
    * as large as the ghost cell width to coarsen.  The source data must have a
    * ghost cell width at least as large as the ghost cell width to coarsen
    * refined to the source (finer) level index space.  Although it is redundant
    * if the coarsen algorithm created the coarsen classes, the routine
    * CoarsenClasses<DIM>::checkCoarsenItem() is also called.
    *
    * If any entries are erroneous an assertion is thrown with a descriptive
    * error message and program halts.
    */
   void initialCheckCoarsenClassItems() const;

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
    * Structures that store coarsen data items.
    */
   tbox::Pointer< xfer::CoarsenClasses<DIM> > d_coarsen_classes;
   int d_number_coarsen_items;
   const typename xfer::CoarsenClasses<DIM>::Data** d_coarsen_items;

   /*!
    * Cached pointers to the coarse, fine, and temporary patch levels.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > d_crse_level;
   tbox::Pointer< hier::PatchLevel<DIM> > d_fine_level;
   tbox::Pointer< hier::PatchLevel<DIM> > d_temp_crse_level;

   /*!
    * Object supporting interface to user-defined spatial data 
    * coarsening operations.
    */
   xfer::CoarsenPatchStrategy<DIM>* d_coarsen_patch_strategy;

   /*!
    * Factory object used to create data transactions when schedule is constructed.
    */
   tbox::Pointer< xfer::CoarsenTransactionFactory<DIM> > d_transaction_factory;

   /*!
    * Cached ratio between source (fine) level and destination (coarse) level.
    */
   hier::IntVector<DIM> d_ratio_between_levels;

   /*!
    * Stored source patch data indices for rapid data allocation/deallocation.
    */
   hier::ComponentSelector d_sources;

   /*!
    * Level-to-level communication schedule between the temporary coarse level
    * and (actual) destination level. 
    */
   tbox::Pointer<tbox::Schedule> d_schedule;

   /*!
    * Boolean indicating whether source data on the coarse temporary level must be
    * filled before coarsening operations (see comments for class constructor in
    * header file), and refine algorithm and schedule needed to peform up these 
    * fill operations.
    */
   bool d_fill_coarse_data;
   tbox::Pointer< RefineAlgorithm<DIM> > d_precoarsen_refine_algorithm;
   tbox::Pointer< RefineSchedule<DIM> > d_precoarsen_refine_schedule;

   /*!
    * Timer objects for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_coarsen_data;
   tbox::Pointer<tbox::Timer> t_gen_sched_n_squared;
   tbox::Pointer<tbox::Timer> t_gen_sched_box_graph;
   tbox::Pointer<tbox::Timer> t_gen_sched_box_tree;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenSchedule.C"
#endif
