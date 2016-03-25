//
// File:	MultiblockCoarsenSchedule.h
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 734 $
// Modified:	$Date: 2005-11-18 10:59:59 -0800 (Fri, 18 Nov 2005) $
// Description:	Coarsening schedule for data transfer between AMR levels
//
 
#ifndef included_mblk_MultiblockCoarsenSchedule
#define included_mblk_MultiblockCoarsenSchedule

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_mblk_MultiblockPatchLevel
#include "MultiblockPatchLevel.h"
#endif
#ifndef included_mblk_MultiblockRefineAlgorithm
#include "MultiblockRefineAlgorithm.h"
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
#ifndef included_xfer_CoarsenClasses
#include "CoarsenClasses.h"
#endif
#ifndef included_mblk_MultiblockCoarsenPatchStrategy
#include "MultiblockCoarsenPatchStrategy.h"
#endif

namespace SAMRAI {
   namespace mblk {

/*!
 * @brief Class MultiblockCoarsenSchedule<DIM> encapsulates the AMR
 * communication pattern to coarsen data from a finer level to a coarser level.
 * 
 * Typically, data is coarsened from the interiors of source patch components 
 * on the source patch level into interiors of destination patch components on  
 * the destination level.  If a coarsen operator has a non-zero ghost cell
 * width, then the source ghost cells must be filled before the coarsen
 * schedule is executed.  The communication schedule is executed by calling
 * member function coarsenData().
 * 
 * Each schedule object is typically created by a coarsen algorithm and
 * represents communication dependencies for a particular configuration
 * of the AMR hierarchy.  The communication schedule is only valid for that
 * particular configuration and must be regenerated when the AMR patch
 * hierarchy changes.  However, as long as the patch levels involved in
 * the creation of the schedule remain unchanged, the schedule may be used
 * for multiple communication cycles.  For more information about creating
 * coarsen schedules, see the MultiblockCoarsenAlgorithm<DIM> header file. 
 *
 * @see mblk::MultiblockCoarsenAlgorithm
 * @see xfer::CoarsenAlgorithm
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::CoarsenClasses
 */
 
template<int DIM>
class MultiblockCoarsenSchedule : public tbox::DescribedClass
{
public:
   /*!
    * @brief Constructor for coarsen schedule from fine level to coarse level
    *
    * Constructor to create a coarsen schedule that coarsens data from 
    * source patch data components on the fine level into the destination patch 
    * data components on the coarse level.  In general, this constructor is 
    * called by a MultiblockCoarsenAlgorithm<DIM> object.  For possible
    * variations on data coarsening, see the Multiblock_CoarsenAlgorithm<DIM>
    * class header information.  
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
    *                          constructed by the calling xfer::CoarsenAlgorithm<DIM>
    *                          object.
    * @param multiblock        Multiblock  hierarchy where the operation
    *                          occurs
    * @param coarsen_strategy  Pointer to a coarsen patch strategy object that
    *                          provides user-defined coarsen operations.  This
    *                          pointer may be null, in which case no
    *                          user-defined coarsen operations will be
    *                          performed.
    * @param refine_strategy   Pointer to a refine patch strategy object that
    *                          provides user-defined coarsen operations.  This
    *                          is needed in specific cases where a user-defined
    *                          coarsen operation requires data to be filled
    *                          on the coarse scratch level prior to execution
    *                          of the coarsening operator.  If such
    *                          functionality is not needed, set this pointer
    *                          to null.
    * @param fill_coarse_data  Boolean indicating whether coarse data should
    *                          be filled before coarsening operations are done.
    * 
    * When assertion checking is active, unrecoverable assertions will result
    * if either patch level pointer, or the refine classes pointer, is null.
    */
   MultiblockCoarsenSchedule(
      tbox::Pointer< MultiblockPatchLevel<DIM> > crse_level,
      tbox::Pointer< MultiblockPatchLevel<DIM> > fine_level,
      const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes,
      tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
      MultiblockCoarsenPatchStrategy<DIM>* coarsen_strategy,
      MultiblockRefinePatchStrategy<DIM>* refine_strategy,
      bool fill_coarse_data);

   /*!
    * @brief The virtual destructor for the schedule releases all internal
    *        storage.
    */
   virtual ~MultiblockCoarsenSchedule<DIM>();

   /*!
    * @brief Execute the stored communication schedule and perform the data
    *        movement.
    */
   void coarsenData() const;

   /*!
    * @brief Return const reference to the pointer to coarsen equivalence
    *        classes used in schedule.
    */
   const tbox::Pointer< xfer::CoarsenClasses<DIM> >& getEquivalenceClasses() const;

private:
   // These two functions are not implemented
   MultiblockCoarsenSchedule(const MultiblockCoarsenSchedule<DIM>&);
   void operator=(const MultiblockCoarsenSchedule<DIM>&);

   /*!
    * @brief set the internal pointer to equivalence classes
    *
    * Utility function to set up local copies of patch data source,
    * destination, etc. indices and necessary data coarsening information
    * stored in the coarsen classes object generated by the coarsen algorithm.
    * An array of coarsen data items is stored locally here to facilitate
    * interaction with transations.
    *
    * @param coarsen_classes   The equivalence classes which store the
    *                           identifiers of the data to be
    */ 
   void setCoarsenItems(
      const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes);

   /*!
    * @brief Utility function to check that patch data has sufficient ghosts.
    *
    * Utility function to check check coarsen items to see whether
    * source and destination patch data components have sufficient ghost
    * cell widths to satisfy the "ghost width to coarsen" functionality
    * described in the xfer::CoarsenAlgorithm<DIM> class header.  Specifically,
    * the destination data must have a ghost cell width at least as large
    * as the ghost cell width to coarsen.  The source data must have a
    * ghost cell width at least as large as the ghost cell width to coarsen
    * refined to the source (finer) level index space.  Although it is
    * redundant if the coarsen algorithm created the coarsen classes, the
    * routine xfer::CoarsenClasses<DIM>::checkCoarsenItem() is also called.
    *
    * If any entries are erroneous an assertion is thrown with a descriptive
    * error message and program halts.
    */
   void initialCheckCoarsenClassItems() const;

   /*!
    * @brief Set up refine algorithm for temporary coarse level filling
    *
    * Set up refine algorithm to fill temporary coarse level before coarsening
    * operations, if necessary.  The associated refine schedule is set in the
    * generateSchedule() routine.
    */
   void setupRefineAlgorithm();

   /*!
    * @brief Generate schedule for moving data from temp storage to destination.
    *
    * Generate communication schedule that moves source patch data from the
    * temporary level into the destination patch data of the destination
    * (coarse) level.
    */
   void generateSchedule();

   /*!
    * @brief Execute the coarsening operation.
    *
    * Coarsen source patch data from the fine patch level into the source patch
    * data on the coarse temporary patch level.
    */
   void coarsenSourceData(MultiblockCoarsenPatchStrategy<DIM> *patch_strategy) const;

   /*!
    * @brief Create a temporary coarse level for data storage during coarsen.
    */
   void generateTemporaryLevel();

   /*
    * Structures that store coarsen data items.
    */
   tbox::Pointer< xfer::CoarsenClasses<DIM> > d_coarsen_classes;
   int d_number_coarsen_items;
   const typename xfer::CoarsenClasses<DIM>::Data** d_coarsen_items;

   /*
    * Cached pointers to the coarse, fine, and temporary patch levels.
    */
   tbox::Pointer< MultiblockPatchLevel<DIM> > d_mblk_crse_level;
   tbox::Pointer< MultiblockPatchLevel<DIM> > d_mblk_fine_level;
   tbox::Pointer< MultiblockPatchLevel<DIM> > d_mblk_temp_crse_level;

   /*
    * Object supporting interface to user-defined spatial data
    * coarsening operations.
    */
   MultiblockCoarsenPatchStrategy<DIM>* d_mblk_coarsen_patch_strategy;

   MultiblockRefinePatchStrategy<DIM>* d_mblk_refine_strategy;

   /*
    * Level-to-level communication schedule between the temporary coarse level
    * and (actual) destination level.
    */
   tbox::Pointer< tbox::Schedule> d_schedule;

   /*
    * Boolean indicating whether source data on the coarse temporary level must
be
    * filled before coarsening operations (see comments for class constructor in    * header file), and refine algorithm and schedule needed to peform up these
    * fill operations.
    */
   bool d_fill_coarse_data;
   tbox::Pointer< MultiblockRefineAlgorithm<DIM> > d_mblk_refine_alg;
   tbox::Pointer< MultiblockRefineSchedule<DIM> > d_mblk_refine_sched;

   hier::ComponentSelector d_sources;

   hier::IntVector<DIM> d_ratio;

   tbox::Pointer< MultiblockPatchHierarchy<DIM> > d_multiblock_hier;
};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockCoarsenSchedule.C"
#endif

#endif
