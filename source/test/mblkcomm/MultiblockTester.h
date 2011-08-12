/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch data communication tests. 
 *
 ************************************************************************/

#ifndef included_MultiblockTester
#define included_MultiblockTester

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"

using namespace std;
using namespace SAMRAI;

class PatchMultiblockTestStrategy;

/**
 * Class MultiblockTester serves as a tool to test data communication operations
 * in SAMRAI, such as coarsening, refining, and filling ghost cells.
 *
 * The functions in this class called from main() are:
 * \begin{enumerate}
 *    - [MultiblockTester(...)] constructor which initializes object state and
 *                            creates patch hierarchy and sets initial data.
 *    - [createRefineSchedule(...)] creates communication schedule for
 *                                      refining data to given level.
 *    - [createCoarsenSchedule(...)] creates communication schedule for
 *                                       coarsening data to given level.
 *    - [performRefineOperations(...)] refines data to given level.
 *    - [performCoarsenOperations(...)] coarsens data to given level.
 * \end{enumerate}
 */

class MultiblockTester:
   public mesh::StandardTagAndInitStrategy,
   public xfer::CoarsenPatchStrategy,
   public xfer::RefinePatchStrategy
{
public:
   /**
    * Constructor performs basic setup operations.
    */
   MultiblockTester(
      const string& object_name,
      const tbox::Dimension& dim,
      tbox::Pointer<tbox::Database>& main_input_db,
      tbox::Pointer<hier::PatchHierarchy>& hierarchy,
      PatchMultiblockTestStrategy* strategy,
      bool do_refine = true,
      bool do_coarsen = false,
      const string& refine_option = "INTERIOR_FROM_SAME_LEVEL");

   /**
    * Destructor is empty.
    */
   ~MultiblockTester();

   /**
    * Return pointer to patch hierarchy on which communication is tested.
    */
   tbox::Pointer<hier::PatchHierarchy> getPatchHierarchy()
   const
   {
      return d_patch_hierarchy;
   }

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void
   registerVariable(
      const tbox::Pointer<hier::Variable> src_variable,
      const tbox::Pointer<hier::Variable> dst_variable,
      const hier::IntVector& src_ghosts,
      const hier::IntVector& dst_ghosts,
      const tbox::Pointer<hier::GridGeometry> xfer_geom,
      const string& operator_name);

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void
   registerVariableForReset(
      const tbox::Pointer<hier::Variable> src_variable,
      const tbox::Pointer<hier::Variable> dst_variable,
      const hier::IntVector& src_ghosts,
      const hier::IntVector& dst_ghosts,
      const tbox::Pointer<hier::GridGeometry> xfer_geom,
      const string& operator_name);

   /**
    * Create communication schedules for refining data to given level.
    */
   void
   createRefineSchedule(
      const int level_number);
   void
   resetRefineSchedule(
      const int level_number);

   /**
    * Create communication schedule for coarsening data to given level.
    */
   void
   createCoarsenSchedule(
      const int level_number);
   void
   resetCoarsenSchedule(
      const int level_number);

   /**
    * Refine data to specified level (or perform interpatch communication
    * on that level).
    */
   void
   performRefineOperations(
      const int level_number);

   /**
    * Coarsen data to specified level.
    */
   void
   performCoarsenOperations(
      const int level_number);

   /**
    * After communication operations are performed, check results.
    */
   bool
   verifyCommunicationResults() const;

   /**
    * Operations needed by GriddingAlgorithm to construct and
    * initialize levels in patch hierarchy.  These operations are
    * pure virtual in GradientDetectorStrategy.
    */
   void
   initializeLevelData(
      const tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const tbox::Pointer<hier::PatchLevel> old_level =
         tbox::Pointer<hier::PatchLevel>(NULL),
      const bool allocate_data = true);

   void
   resetHierarchyConfiguration(
      const tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const int coarsest_level,
      const int finest_level);

   void
   applyGradientDetector(
      const tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * These routines pass off physicial boundary and pre/postprocess
    * coarsen/refine operations to patch data test object.  They are
    * pure virtual in RefinePatchStrategy and CoarsenPatchStrategy.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw);

   /*!
    * Set the ghost data at a multiblock singularity.
    */
   void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      const hier::Connector& dst_to_encon,
      const double fill_time,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box);

   hier::IntVector
   getRefineOpStencilWidth() const;

   void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   hier::IntVector
   getCoarsenOpStencilWidth() const;

   void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   double getLevelDt(
      const tbox::Pointer<hier::PatchLevel> level,
      const double dt_time,
      const bool initial_time)
   {
      (void)level;
      (void)dt_time;
      (void)initial_time;
      return 0.0;
   }

   /*
    * Construct patch hierarchy and initialize data prior to tests.
    */
   void
   setupHierarchy(
      tbox::Pointer<tbox::Database> main_input_db,
      tbox::Pointer<mesh::StandardTagAndInitialize> cell_tagger);

private:
   /*
    * Object name for error reporting.
    */
   string d_object_name;

   const tbox::Dimension d_dim;

   /*
    * Object supplying operatins for particular patch data test.
    */
   PatchMultiblockTestStrategy* d_data_test_strategy;

   /*
    * Booleans to indicate whether refine or coarsen is operation to test.
    */
   bool d_do_refine;
   bool d_do_coarsen;

   /*
    * String name for refine option; ; i.e., source of interior patch
    * data on refined patches.  Options are "INTERIOR_FROM_SAME_LEVEL"
    * and "INTERIOR_FROM_COARSER_LEVEL".
    */
   string d_refine_option;

   /*
    * Patch hierarchy on which tests occur.
    */
   tbox::Pointer<hier::PatchHierarchy> d_patch_hierarchy;

   /*
    * Dummy time stamp for all data operations.
    */
   double d_fake_time;

   /*
    * The MultiblockTester uses two variable contexts for each variable.
    * The "source", and "destination" contexts indicate the source
    * and destination patch data for the transfer operation.
    *
    * The "refine_scratch" context is used for managing scratch
    * space during refine operations.
    */
   tbox::Pointer<hier::VariableContext> d_source;
   tbox::Pointer<hier::VariableContext> d_destination;
   tbox::Pointer<hier::VariableContext> d_refine_scratch;

   tbox::Pointer<hier::VariableContext> d_reset_source;
   tbox::Pointer<hier::VariableContext> d_reset_destination;
   tbox::Pointer<hier::VariableContext> d_reset_refine_scratch;

   /*
    * Component selector for allocation/deallocation of variable data.
    */
   hier::ComponentSelector d_patch_data_components;

   /*
    * Refine/Coarsen algorithm and schedules for testing communication
    * among levels in the patch hierarchy.
    */

   tbox::Pointer<xfer::RefineAlgorithm> d_refine_algorithm;
   tbox::Pointer<xfer::CoarsenAlgorithm> d_coarsen_algorithm;

   xfer::RefineAlgorithm d_reset_refine_algorithm;
   xfer::CoarsenAlgorithm d_reset_coarsen_algorithm;

   tbox::Pointer<xfer::RefineAlgorithm> d_mblk_refine_alg;

   bool d_is_reset;

   tbox::Array<tbox::Pointer<xfer::RefineSchedule> >
   d_refine_schedule;
   tbox::Array<tbox::Pointer<xfer::CoarsenSchedule> >
   d_coarsen_schedule;

};

#endif
