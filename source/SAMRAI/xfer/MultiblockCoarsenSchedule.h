/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Coarsening schedule for data transfer between AMR levels 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockCoarsenSchedule
#define included_xfer_MultiblockCoarsenSchedule

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/xfer/CoarsenClasses.h"
#include "SAMRAI/xfer/CoarsenTransactionFactory.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/MultiblockCoarsenPatchStrategy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class MultiblockCoarsenSchedule encapsulates the AMR
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
 * coarsen schedules, see the MultiblockCoarsenAlgorithm header file.
 *
 * @see xfer::MultiblockCoarsenAlgorithm
 * @see xfer::CoarsenAlgorithm
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::CoarsenClasses
 */

class MultiblockCoarsenSchedule:public tbox::DescribedClass
{
public:
   /*!
    * Static function to set box intersection algorithm to use during
    * schedule construction for all CoarsenSchedule objects.
    * If this method is not called, the default will be used.
    *
    * @param  method   string identifying box intersection method.  Valid
    *                  choices are:  "DLBG" (default case),
    *                  and "ORIG_NSQUARED".   More details can be found below
    *                  in the comments for the generateSchedule() routine.
    *
    * If an invalid string is passed, an unrecoverable error will result.
    * The "ORIG_NSQUARED" option only works for problems with one block.
    */
   static void
   setScheduleGenerationMethod(
      const std::string& method);

   /*!
    * @brief Constructor for coarsen schedule from fine level to coarse level
    *
    * Constructor to create a coarsen schedule that coarsens data from
    * source patch data components on the fine level into the destination patch
    * data components on the coarse level.  In general, this constructor is
    * called by a MultiblockCoarsenAlgorithm object.  For possible
    * variations on data coarsening, see the Multiblock_CoarsenAlgorithm
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
    *                          constructed by the calling xfer::CoarsenAlgorithm
    *                          object.
    * @param hierarchy         Multiblock  hierarchy where the operation
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
   explicit MultiblockCoarsenSchedule(
      tbox::Pointer<hier::PatchLevel> crse_level,
      tbox::Pointer<hier::PatchLevel> fine_level,
      const tbox::Pointer<xfer::CoarsenClasses> coarsen_classes,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      tbox::Pointer<xfer::CoarsenTransactionFactory> transaction_factory,
      MultiblockCoarsenPatchStrategy* coarsen_strategy,
      RefinePatchStrategy* refine_strategy,
      bool fill_coarse_data);

   /*!
    * @brief The virtual destructor for the schedule releases all internal
    *        storage.
    */
   virtual ~MultiblockCoarsenSchedule();

   /*!
    * @brief Execute the stored communication schedule and perform the data
    *        movement.
    */
   void
   coarsenData() const;

   /*!
    * @brief Return const reference to the pointer to coarsen equivalence
    *        classes used in schedule.
    */
   const tbox::Pointer<xfer::CoarsenClasses>&
   getEquivalenceClasses() const;

private:
   // These two functions are not implemented
   MultiblockCoarsenSchedule(
      const MultiblockCoarsenSchedule&);
   void
   operator = (
      const MultiblockCoarsenSchedule&);

   //! @brief Shorthand typedef.
   typedef hier::LocalId LocalId;
   //! @brief Shorthand typedef.
   typedef hier::Connector::NeighborSet NeighborSet;
   //! @brief Mapping from a (potentially remote) MappedBox to a set of neighbors.
   typedef std::map<hier::MappedBox, NeighborSet> FullNeighborhoodSet;

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
   void
   setCoarsenItems(
      const tbox::Pointer<xfer::CoarsenClasses> coarsen_classes);

   /*!
    * @brief Utility function to clear local copies of coarsen items.
    */
   void
   clearCoarsenItems();

   /*!
    * @brief Utility function to check that patch data has sufficient ghosts.
    *
    * Utility function to check check coarsen items to see whether
    * source and destination patch data components have sufficient ghost
    * cell widths to satisfy the "ghost width to coarsen" functionality
    * described in the xfer::CoarsenAlgorithm class header.  Specifically,
    * the destination data must have a ghost cell width at least as large
    * as the ghost cell width to coarsen.  The source data must have a
    * ghost cell width at least as large as the ghost cell width to coarsen
    * refined to the source (finer) level index space.  Although it is
    * redundant if the coarsen algorithm created the coarsen classes, the
    * routine xfer::CoarsenClasses::checkCoarsenItem() is also called.
    *
    * If any entries are erroneous an assertion is thrown with a descriptive
    * error message and program halts.
    */
   void
   initialCheckCoarsenClassItems() const;

   /*!
    * @brief Set up refine algorithm for temporary coarse level filling
    *
    * Set up refine algorithm to fill temporary coarse level before coarsening
    * operations, if necessary.  The associated refine schedule is set in the
    * generateSchedule() routine.
    */
   void
   setupRefineAlgorithm();

   /*!
    * @brief Generate schedule for moving data from temp storage to destination.
    *
    * Generate communication schedule that moves source patch data from the
    * temporary level into the destination patch data of the destination
    * (coarse) level.
    */
   void
   generateSchedule();

   /*!
    * @brief This version of the schedule generation procedure
    * uses N^2 algorithms to determine box intersections;
    * i.e., the original SAMRAI implementation which checks every
    * box against every other.
    */
   void
   generateScheduleNSquared();

   /*!
    * @brief This version of the schedule generation procedure uses the DLBG
    * data to determine which source patches contribute
    * data to each destination patch and to compute unfilled_boxes.
    */
   void
   generateScheduleDLBG();

   /*!
    * @brief Execute the coarsening operation.
    *
    * Coarsen source patch data from the fine patch level into the source patch
    * data on the coarse temporary patch level.
    */
   void
   coarsenSourceData(
      MultiblockCoarsenPatchStrategy* patch_strategy) const;

   /*!
    * @brief Calculate the maximum ghost cell width to grow boxes to check
    * for overlaps.
    */
   hier::IntVector
   getMaxGhostsToGrow() const;

   /*!
    * @brief Function that constructs schedule transactions that
    * move data from source patch on source level to destination patch
    * on destination level.
    */
   void
   constructScheduleTransactions(
      tbox::Pointer<hier::PatchLevel> dst_level,
      const hier::MappedBox& dst_mapped_box,
      tbox::Pointer<hier::PatchLevel> src_level,
      const hier::MappedBox& src_mapped_box);

   void
   constructScheduleTransactions(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level);

   /*!
    * @brief Constructs the RefineSchedule that is used
    * to move coarsened data from the temporary coarse level to the
    * destination level.
    */
   void
   constructDestinationLevelFillSchedule(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level);

   /*!
    * @brief Restructure a src_to_dst Connector to arange edges by dst
    * major order, and apply shifts to cancel out any periodic shifts
    * on dst mapped_boxes.
    */
   void
   restructureNeighborhoodSetsByDstNodes(
      const hier::Connector& cnect,
      FullNeighborhoodSet& full_inverted_edges) const;

   /*!
    * @brief Create a temporary coarse level for data storage during coarsen.
    */
   void
   generateTemporaryLevel();

   /*
    * @brief Get the overlap Connector registered for the given base
    * and head.
    */
   const hier::Connector *
   getOverlapConnector(
      const hier::MappedBoxLevel& base,
      const hier::MappedBoxLevel& head,
      const hier::IntVector& min_gcw) const;

   /*
    * @brief Similar to getOverlapConnector() but requires that the
    * Connector is found.
    */
   const hier::Connector *
   getOverlapConnector_strict(
      const hier::MappedBoxLevel& base,
      const hier::MappedBoxLevel& head,
      const hier::IntVector& min_gcw) const;

   /*!
    * Selects algorithm used to generate communication schedule.
    */
   static std::string s_schedule_generation_method;

   /*
    * Structures that store coarsen data items.
    */
   tbox::Pointer<xfer::CoarsenClasses> d_coarsen_classes;
   int d_number_coarsen_items;
   const xfer::CoarsenClasses::Data** d_coarsen_items;

   /*
    * Cached pointers to the coarse, fine, and temporary patch levels.
    */
   tbox::Pointer<hier::PatchLevel> d_mblk_crse_level;
   tbox::Pointer<hier::PatchLevel> d_mblk_fine_level;
   tbox::Pointer<hier::PatchLevel> d_mblk_temp_crse_level;

   /*!
    * @brief Fine-to-coarse connector, NULL if data not provided.
    */
   const hier::Connector** d_fine_to_coarse;

   /*!
    * @brief Coarse-to-fine connector, NULL if data not provided.
    */
   const hier::Connector** d_coarse_to_fine;
   /*!
    * @brief Connector from temporary (coarsened fine) mapped_box_level
    * to coarse mapped_box_level.
    */
   tbox::Array<hier::Connector> d_temp_to_coarse;
   /*!
    * @brief Connector from coarse mapped_box_level to temporary
    * (coarsened fine) mapped_box_level.
    */
   tbox::Array<hier::Connector> d_coarse_to_temp;

   /*
    * Object supporting interface to user-defined spatial data
    * coarsening operations.
    */
   MultiblockCoarsenPatchStrategy* d_mblk_coarsen_patch_strategy;

   RefinePatchStrategy* d_mblk_refine_strategy;

   /*!
    * Factory object used to create data transactions when schedule is constructed.
    */
   tbox::Pointer<xfer::CoarsenTransactionFactory> d_transaction_factory;

   /*
    * Level-to-level communication schedule between the temporary coarse level
    * and (actual) destination level.
    */
   tbox::Pointer<tbox::Schedule> d_schedule;

   /*
    * Boolean indicating whether source data on the coarse temporary level must
    * be
    * filled before coarsening operations (see comments for class constructor in    * header file), and refine algorithm and schedule needed to peform up these
    * fill operations.
    */
   bool d_fill_coarse_data;
   tbox::Pointer<RefineAlgorithm> d_mblk_fill_coarse_data_alg;
   tbox::Pointer<RefineSchedule> d_mblk_fill_coarse_data_sched;

   tbox::Pointer<RefineAlgorithm> d_mblk_fill_dst_alg;
   tbox::Pointer<RefineSchedule> d_mblk_fill_dst_sched;

   hier::ComponentSelector d_sources;

   hier::IntVector d_ratio_between_levels;

   tbox::Pointer<hier::PatchHierarchy> d_mblk_hierarchy;

   /*!
    * Timer objects for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_coarsen_data;
   tbox::Pointer<tbox::Timer> t_gen_sched_n_squared;
   tbox::Pointer<tbox::Timer> t_gen_sched_dlbg;
   tbox::Pointer<tbox::Timer> t_modify_connector;
   tbox::Pointer<tbox::Timer> t_invert_edges;
};

}
}

#endif
