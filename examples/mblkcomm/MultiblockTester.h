// 
// File:        MultiblockTester.h
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 1.7 $
// Modified:    $Date: 2004/03/08 21:12:08 $
// Description: Manager class for patch data communication tests.
//

#ifndef included_MultiblockTester
#define included_MultiblockTester

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_BasePatchHierarchy
#include "BasePatchHierarchy.h"
#endif
#ifndef included_hier_BasePatchLevel
#include "BasePatchLevel.h"
#endif
#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_xfer_CoarsenAlgorithm
#include "CoarsenAlgorithm.h"
#endif
#ifndef included_xfer_CoarsenPatchStrategy
#include "CoarsenPatchStrategy.h"
#endif
#ifndef included_xfer_CoarsenSchedule
#include "CoarsenSchedule.h"
#endif
#ifndef included_tbox_Database
#include "Database.h"
#endif
#ifndef included_xfer_Geometry
#include "Geometry.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_mblk_MultiblockCoarsenAlgorithm
#include "MultiblockCoarsenAlgorithm.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_mblk_MultiblockRefinePatchStrategy
#include "MultiblockRefinePatchStrategy.h"
#endif
#ifndef included_mblk_MultiblockRefineAlgorithm
#include "MultiblockRefineAlgorithm.h"
#endif
#ifndef included_mblk_MultiblockRefineSchedule
#include "MultiblockRefineSchedule.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "Pointer.h"
#endif
#ifndef included_xfer_RefineAlgorithm
#include "RefineAlgorithm.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif
#ifndef included_xfer_RefineSchedule
#include "RefineSchedule.h"
#endif
#ifndef included_mesh_StandardTagAndInitialize
#include "StandardTagAndInitialize.h"
#endif
#ifndef included_mesh_StandardTagAndInitStrategy
#include "StandardTagAndInitStrategy.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

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

class MultiblockTester :
   public mesh::StandardTagAndInitStrategy<NDIM>,
   public xfer::CoarsenPatchStrategy<NDIM>,
   public mblk::MultiblockRefinePatchStrategy<NDIM>
{
public:

   /**
    * Constructor performs basic setup operations.
    */
   MultiblockTester(const string& object_name,
              tbox::Pointer<tbox::Database> main_input_db,
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
   tbox::Pointer< mblk::MultiblockPatchHierarchy<NDIM> > getPatchHierarchy()
   const
   {
      return(d_patch_hierarchy);
   }

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void registerVariable(
      const tbox::Pointer< hier::Variable<NDIM> > src_variable,
      const tbox::Pointer< hier::Variable<NDIM> > dst_variable,
      const hier::IntVector<NDIM>& src_ghosts,
      const hier::IntVector<NDIM>& dst_ghosts,
      const tbox::Pointer< xfer::Geometry<NDIM> > xfer_geom,
      const string& operator_name);

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void registerVariableForReset(
      const tbox::Pointer< hier::Variable<NDIM> > src_variable,
      const tbox::Pointer< hier::Variable<NDIM> > dst_variable,
      const hier::IntVector<NDIM>& src_ghosts,
      const hier::IntVector<NDIM>& dst_ghosts,
      const tbox::Pointer< xfer::Geometry<NDIM> > xfer_geom,
      const string& operator_name);

   /**
    * Create communication schedules for refining data to given level.
    */
   void createRefineSchedule(const int level_number);
   void resetRefineSchedule(const int level_number);

   /**
    * Create communication schedule for coarsening data to given level.
    */
   void createCoarsenSchedule(const int level_number);
   void resetCoarsenSchedule(const int level_number);

   /**
    * Refine data to specified level (or perform interpatch communication
    * on that level).
    */
   void performRefineOperations(const int level_number);

   /**
    * Coarsen data to specified level.
    */
   void performCoarsenOperations(const int level_number);

   /**
    * After communication operations are performed, check results.
    */
   void verifyCommunicationResults() const;

   /**
    * Operations needed by GriddingAlgorithm to construct and
    * initialize levels in patch hierarchy.  These operations are
    * pure virtual in GradientDetectorStrategy.
    */

   void initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
      const int level_number,
      const double init_time,
      const bool can_be_refined,
      const bool initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level =
         tbox::Pointer< hier::BasePatchLevel<NDIM> >(NULL),
      const bool allocate_data = true);

   void resetHierarchyConfiguration(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   void applyGradientDetector(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
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
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const hier::IntVector<NDIM>& gcw);

   /*!
    * Set the ghost data at a multiblock singularity.
    */
   void fillSingularityBoundaryConditions(
      hier::Patch<NDIM>& patch,
      tbox::List<mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
         singularity_patches,
      const double fill_time,
      const hier::Box<NDIM>& fill_box,
      const hier::BoundaryBox<NDIM>& boundary_box);

   hier::IntVector<NDIM> getRefineOpStencilWidth() const;

   void preprocessRefine(hier::Patch<NDIM>& fine,
                         const hier::Patch<NDIM>& coarse,
                         const hier::Box<NDIM>& fine_box,
                         const hier::IntVector<NDIM>& ratio);

   void postprocessRefine(hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio);

   hier::IntVector<NDIM> getCoarsenOpStencilWidth() const;

   void preprocessCoarsen(hier::Patch<NDIM>& coarse,
                          const hier::Patch<NDIM>& fine,
                          const hier::Box<NDIM>& coarse_box,
                          const hier::IntVector<NDIM>& ratio);

   void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                           const hier::Patch<NDIM>& fine,
                           const hier::Box<NDIM>& coarse_box,
                           const hier::IntVector<NDIM>& ratio);

   double getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                     const double dt_time,
                     const bool initial_time)
   {
       (void) level;
       (void) dt_time;
       (void) initial_time;
       return (0.0);
   }

   /*
    * Construct patch hierarchy and initialize data prior to tests.
    */
   void setupHierarchy(
      tbox::Pointer<tbox::Database> main_input_db, 
      tbox::Pointer< mesh::StandardTagAndInitialize<NDIM> > cell_tagger);

private:

   /*
    * Object name for error reporting.
    */
   string d_object_name;

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
   tbox::Pointer< mblk::MultiblockPatchHierarchy<NDIM> > d_patch_hierarchy;

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

   tbox::Pointer< xfer::RefineAlgorithm<NDIM> >  d_refine_algorithm;
   tbox::Pointer< xfer::CoarsenAlgorithm<NDIM> > d_coarsen_algorithm;

   xfer::RefineAlgorithm<NDIM>  d_reset_refine_algorithm;
   xfer::CoarsenAlgorithm<NDIM> d_reset_coarsen_algorithm;

   tbox::Pointer <mblk::MultiblockRefineAlgorithm<NDIM> > d_mblk_refine_alg;

   bool d_is_reset;

   tbox::Array< tbox::Pointer< mblk::MultiblockRefineSchedule<NDIM> > >
      d_refine_schedule;
   tbox::Array< tbox::Pointer< xfer::CoarsenSchedule<NDIM> > >
      d_coarsen_schedule;


};

#endif
