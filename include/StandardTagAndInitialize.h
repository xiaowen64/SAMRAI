//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/gridding/StandardTagAndInitialize.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Gridding routines and params for Richardson Extrapolation.
//

#ifndef included_mesh_StandardTagAndInitialize
#define included_mesh_StandardTagAndInitialize

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_mesh_StandardTagAndInitStrategy
#include "StandardTagAndInitStrategy.h"
#endif
#ifndef included_mesh_TagAndInitializeStrategy
#include "TagAndInitializeStrategy.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


namespace SAMRAI {
    namespace mesh {

/*!
 * Class StandardTagAndInitialize<DIM> defines an implementation
 * for level initialization and cell tagging routines needed by
 * the GriddingAlgorithm<DIM> class.  This class is derived from
 * the abstract base class TagAndInitializeStrategy<DIM>. It invokes 
 * problem-specific level initialization routines after AMR patch 
 * hierarchy levels change and routines for tagging cells for refinement
 * using one (or more) of the following methods:
 *
 *   - Gradient Detection
 *   - Richardson Extrapolation
 *   - Explicitly defined refine boxes
 *
 * It is possible to use combinations of these three methods (e.g., 
 * use gradient detection, Richardson extrapolation, and static refine boxes
 * at the same time). The order in which they are executed is fixed (
 * Richardson extrapolation first, gradient detection second, and refine
 * boxes third).  An input entry for this class is optional.
 * If none is provided, the class will by default not use any criteria  
 * to tag cells for refinement.
 * 
 * Required input keys and data types: NONE
 *
 * Optional input keys, data types, and defaults:
 * 
 *
 *    - \b    tagging_method   
 *       string array specification of the type of cell-tagging used.  Valid
 *       choices include:
 *          - ``GRADIENT_DETECTOR''
 *          - ``RICHARDSON_EXTRAPOLATION''
 *          - ``REFINE_BOXES'' 
 *          - ``RICHARDSON_EXTRAPOLATION'', ``GRADIENT_DETECTOR'', 
 *                ``REFINE_BOXES''
 *       (i.e. a combination of any or all of the above - the choices may 
 *       be placed in any order). If no input is given, no tagging will be 
 *       performed.
 *
 *    - \b    input section describing the refine boxes for each level.
 *      (@see mesh::TagAndInitializeStrategy for details on format)
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *    tagging_method = "GRADIENT_DETECTOR", "REFINE_BOXES"
 *    <refine boxes input> (@see mesh::TagAndInitializeStrategy)
 *    
 * \endverbatim
 *
 * This class supplies the routines for tagging cells 
 * and invoking problem-specific level initialization routines after AMR 
 * patch hierarchy levels change.  A number of problem-specific operations 
 * are performed in the StandardTagAndInitStrategy<DIM> 
 * data member, for which methods are specified in a derived subclass.
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::GriddingAlgorithm
 * @see mesh::StandardTagAndInitStrategy
 */

template<int DIM> class StandardTagAndInitialize 
: 
public TagAndInitializeStrategy<DIM>  
{
public:
   /*!
    * Constructor for StandardTagAndInitialize<DIM> which
    * may read inputs from the provided input_db.  If no input 
    * database is provided, the class interprets that no tagging
    * is desired so no cell-tagging will be performed.   
    */
   StandardTagAndInitialize(
      const std::string& object_name,
      StandardTagAndInitStrategy<DIM>* tag_strategy,
      tbox::Pointer<tbox::Database> input_db = tbox::Pointer<tbox::Database>(NULL));

   /*!
    * Virtual destructor for StandardTagAndInitialize<DIM>.
    */
   virtual ~StandardTagAndInitialize<DIM>();

   /*!
    * Specifies whether the chosen method advances the solution data 
    * in the regridding process (Richardson extrapolation does, the
    * others will not).
    */
   bool usesTimeIntegration() const;

   /*!
    * Return coarsen ratio used for applying cell tagging. An error
    * coarsen ratio other than 2 or 3 will throw an error.
    */
   int getErrorCoarsenRatio() const;

   /*!
    * Some restrictions may be placed on the coarsen ratio used for
    * cell tagging.  Check these here.
    */
   void
   checkCoarsenRatios(const tbox::Array< hier::IntVector<DIM> >& ratio_to_coarser);

   /*!
    * Pass the request to initialize the data on a new level in the 
    * hierarchy to the StandardTagAndInitStrategy<DIM> data member. Required
    * arguments specify the grid hierarchy, level number being initialized,
    * simulation time at which the data is initialized, whether the level
    * can be refined, and whether it is the initial time.  Optional arguments
    * include an old level, from which data may be used to initialize this
    * level, and a flag that indicates whether data on the initialized level 
    * must first be allocated.  For more information on the operations that 
    * must be performed, see the 
    * TagAndInitializeStrategy<DIM>::initializeLevelData() method.
    */
   void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                       const int level_number,
                       const double init_data_time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<DIM> > old_level = 
                             tbox::Pointer< hier::BasePatchLevel<DIM> >(NULL),
		       const bool allocate_data = true);

   /*!
    * Pass the request to reset information that depends on the hierarchy
    * configuration to the StandardTagAndInitStrategy<DIM> data member.  
    * For more information on the operations that must be performed, see 
    * the TagAndInitializeStrategy<DIM>::resetHierarchyConfiguration()
    * method.
    */
   void resetHierarchyConfiguration(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    * Certain cases may require pre-processing of error estimation data
    * before tagging cells, which is handled by this method.  For more 
    * information on the operations that must be performed, see the
    * TagAndInitializeStrategy<DIM>::preprocessErrorEstimation()
    * method
    */
   void preprocessErrorEstimation(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time, 
      const double regrid_start_time,
      const bool initial_time);

   /*!
    * Pass the request to set tags on the given level where refinement of 
    * that level should occur.  Gradient detection, Richardson extrapolation,
    * and tagging on static refine boxes is performed here. 
    *
    * For more information on the operations that must be performed, see the 
    * TagAndInitializeStrategy<DIM>::tagCellsForRefinement() routine.
    */
   void tagCellsForRefinement(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > level,
      const int level_number,
      const double regrid_time,
      const int tag_index,
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined,
      const double regrid_start_time = 0);
			     
   /*!
    * Return true if boxes for coarsest hierarchy level are not appropriate
    * for gridding strategy.  Otherwise, return false.  If false is returned,
    * it is useful to provide a detailed explanatory message describing the
    * problems with the boxes.
    */
   bool coarsestLevelBoxesOK(const hier::BoxArray<DIM>& boxes) const;

   /*!
    * Return whether refinement is being performed using ONLY 
    * user-supplied refine boxes.  If any method is used that invokes
    * tagging, this will return false.
    */
   bool refineUserBoxInputOnly() const;

   /*!
    * Turn on gradient detector to tag cells for refinement. 
    */
   void turnOnGradientDetector();

   /*!
    * Turn off gradient detector.
    */
   void turnOffGradientDetector();

   /*!
    * Turn on Richardson extrapolation to tag cells for refinement. 
    */
   void turnOnRichardsonExtrapolation();

   /*!
    * Turn off Richardson extrapolation.
    */
   void turnOffRichardsonExtrapolation();

   /*!
    * Turn on static refine box regions where refinement should occur. 
    */
   void turnOnRefineBoxes();

   /*!
    * Turn off static refine box regions. 
    */
   void turnOffRefineBoxes();

   /*!
    * Read input values, indicated above, from given database. 
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

private:
   /*
    * Apply preprocessing for Richardson extrapolation.
    */
   void preprocessRichardsonExtrapolation(
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time, 
      const double regrid_start_time,
      const bool initial_time);

   /*
    * Apply Richardson extrapolation algorithm.
    */
   void tagCellsUsingRichardsonExtrapolation(
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      const double regrid_start_time,
      const int tag_index, 
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined);

   /*
    * Object name.
    */
   std::string d_object_name;

   /*
    * Booleans specifying the tagging method.  Any combination of the
    * three methods may be used.
    */
   bool d_use_refine_boxes;
   bool d_use_gradient_detector;
   bool d_use_richardson_extrapolation;

   /*
    * Concrete object that supplies problem-specific initialization
    * and regridding operations.
    */
   StandardTagAndInitStrategy<DIM>* d_tag_strategy;

   /*
    * The error_coarsen_ratio used for all levels in the hierarchy.
    * If Richardson extrapolation is not used, the error coarsen ratio
    * is 1.  If Richardson extrapolation is used, the error coarsen ratio
    * is set in the method coarsestLevelBoxesOK().
    */
   int d_error_coarsen_ratio;

   /*
    * tbox::Array of patch levels containing coarsened versions of the patch
    * levels, for use with Richardson extrapolation.
    */
   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > > d_rich_extrap_coarsened_levels;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "StandardTagAndInitialize.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "StandardTagAndInitialize.C"
#endif
