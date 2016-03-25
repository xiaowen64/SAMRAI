// 
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/communication/CommTester.h $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Manager class for patch data communication tests.
//

#ifndef included_CommTester
#define included_CommTester

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
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
#include "tbox/Database.h"
#endif
#ifndef included_xfer_Geometry
#include "Geometry.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchDataTestStrategy
#include "PatchDataTestStrategy.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
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
#ifndef included_String
#include <string>
using namespace std;
#define included_String
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

namespace SAMRAI {

/**
 * Class CommTester serves as a tool to test data communication operations 
 * in SAMRAI, such as coarsening, refining, and filling ghost cells.  
 * 
 * The functions in this class called from main() are:
 * \begin{enumerate}
 *    - [CommTester(...)] constructor which initializes object state and 
 *                            creates patch hierarchy and sets initial data.
 *    - [createRefineSchedule(...)] creates communication schedule for
 *                                      refining data to given level.
 *    - [createCoarsenSchedule(...)] creates communication schedule for
 *                                       coarsening data to given level.
 *    - [performRefineOperations(...)] refines data to given level.
 *    - [performCoarsenOperations(...)] coarsens data to given level.
 * \end{enumerate}
 */

class CommTester :
   public mesh::StandardTagAndInitStrategy<NDIM>,
   public xfer::RefinePatchStrategy<NDIM>,
   public xfer::CoarsenPatchStrategy<NDIM>
{
public:

   /**
    * Constructor performs basic setup operations.
    */
   CommTester(const string& object_name,
              tbox::Pointer<tbox::Database> main_input_db,
              PatchDataTestStrategy* strategy,
              bool do_refine = true,
              bool do_coarsen = false,
              const string& refine_option = "INTERIOR_FROM_SAME_LEVEL");

   /**
    * Destructor is empty.
    */
   ~CommTester();

   /**
    * Return pointer to patch hierarchy on which communication is tested.
    */
   tbox::Pointer<hier::PatchHierarchy<NDIM> > getPatchHierarchy() const
   {
      return(d_patch_hierarchy);
   }

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void registerVariable(const tbox::Pointer<hier::Variable<NDIM> > src_variable,
                         const tbox::Pointer<hier::Variable<NDIM> > dst_variable,
                         const hier::IntVector<NDIM>& src_ghosts,
                         const hier::IntVector<NDIM>& dst_ghosts,
                         const tbox::Pointer<xfer::Geometry<NDIM> > xfer_geom,
                         const string& operator_name);

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void registerVariableForReset(const tbox::Pointer<hier::Variable<NDIM> > src_variable,
                                 const tbox::Pointer<hier::Variable<NDIM> > dst_variable,
                                 const hier::IntVector<NDIM>& src_ghosts,
                                 const hier::IntVector<NDIM>& dst_ghosts,
                                 const tbox::Pointer<xfer::Geometry<NDIM> > xfer_geom,
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
   bool verifyCommunicationResults() const;

   /**
    * Operations needed by mesh::GriddingAlgorithm<NDIM> to construct and
    * initialize levels in patch hierarchy.  These operations are
    * pure virtual in GradientDetectorStrategy.
    */

   void initializeLevelData(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                            const int level_number,
                            const double init_time,
                            const bool can_be_refined,
                            const bool initial_time,
                            const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level =
                                  tbox::Pointer<hier::BasePatchLevel<NDIM> >(NULL),
                            const bool allocate_data = true);

   void resetHierarchyConfiguration(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                                    const int coarsest_level,
                                    const int finest_level);
   /**
    * These routines pass off physicial boundary and pre/postprocess 
    * coarsen/refine operations to patch data test object.  They are
    * pure virtual in RefinePatchStrategy and CoarsenPatchStrategy.
    */
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const hier::IntVector<NDIM>& gcw);

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

   double getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
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
   void setupHierarchy(tbox::Pointer<tbox::Database> main_input_db, 
                       tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > cell_tagger);

private:
   /*
    * Object name for error reporting.
    */
   string d_object_name;

   /*
    * Object supplying operatins for particular patch data test.
    */
   PatchDataTestStrategy* d_data_test_strategy;

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
    *hier::Patch<NDIM> hierarchy on which tests occur.
    */
   tbox::Pointer<hier::PatchHierarchy<NDIM> > d_patch_hierarchy;

   /*
    * Dummy time stamp for all data operations.
    */
   double d_fake_time;

   /*
    * The CommTester uses two variable contexts for each variable.
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

   xfer::RefineAlgorithm<NDIM>  d_refine_algorithm;
   xfer::CoarsenAlgorithm<NDIM> d_coarsen_algorithm;

   xfer::RefineAlgorithm<NDIM>  d_reset_refine_algorithm;
   xfer::CoarsenAlgorithm<NDIM> d_reset_coarsen_algorithm;

   bool d_is_reset;

   tbox::Array< tbox::Pointer<xfer::RefineSchedule<NDIM> > > d_refine_schedule;
   tbox::Array< tbox::Pointer<xfer::CoarsenSchedule<NDIM> > > d_coarsen_schedule;

};

}
#endif
