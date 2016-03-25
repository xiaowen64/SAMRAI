//
// File:        MethodOfLinesPatchStrategy.h
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to application-specific patch functions in support
//              Method of Lines integration algorithm   
//

#ifndef included_algs_MethodOfLinesPatchStrategy
#define included_algs_MethodOfLinesPatchStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif
#ifndef included_xfer_CoarsenPatchStrategy
#include "CoarsenPatchStrategy.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace algs {

template<int DIM> class MethodOfLinesIntegrator;

/**
 * Class MethodOfLinesPatchStrategy<DIM> is an abstract type defining the 
 * interface for operations invoked during the integration routines defined 
 * in the MethodOfLinesIntegrator<DIM> class.  This class is derived from 
 * the xfer::RefinePatchStrategy<DIM> and xfer::CoarsenPatchStrategy<DIM> abstract 
 * base classes.  These base classes define the interfaces for user-defined 
 * interlevel data refining and coarsening operations and the specification 
 * of physical boundary conditions.
 *
 * @see algs::MethodOfLinesIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

template<int DIM> class MethodOfLinesPatchStrategy 
:
   public xfer::RefinePatchStrategy<DIM>,
   public xfer::CoarsenPatchStrategy<DIM> 
{
public:
   /*!
    * Blank constructor for MethodOfLinesPatchStrategy<DIM>.
    */
   MethodOfLinesPatchStrategy();

   /*!
    * Virtual destructor for MethodOfLinesPatchStrategy<DIM>.
    */
   virtual ~MethodOfLinesPatchStrategy<DIM>() = 0;

   /*!
    * Register variables specific to the problem to be solved with the 
    * integrator using the registerVariable function.  This 
    * defines the way data for each quantity will be manipulated on the 
    * patches.  For more information, refer to   
    * MethodOfLinesIntegrator<DIM>::registerVariable.
    */
   virtual void 
   registerModelVariables(MethodOfLinesIntegrator<DIM>* integrator) = 0; 

   /*!
    * Set the initial data on a patch interior (i.e., NO GHOST CELLS).  
    * Setting "initial_time" true will initialize data at time. Setting it
    * false will interpolate data from the appropriate coarser patch. 
    */
   virtual void initializeDataOnPatch(hier::Patch<DIM>& patch,
                                      const double time,
                                      const bool initial_time) const = 0;

   /*!
    * Compute the stable time increment for a patch.   
    */
   virtual double computeStableDtOnPatch(hier::Patch<DIM>& patch,
               const double time) const = 0;

   /*!
    * Advance a single Runge Kutta step.
    * 
    * @param patch patch that RK step is being applied
    * @param dt    timestep
    * @param alpha_1 first coefficient applied in the RK step
    * @param alpha_2 second coefficient
    * @param beta    third coefficient
    */
   virtual void singleStep(hier::Patch<DIM>& patch,
                           const double dt, 
                           const double alpha_1,
                           const double alpha_2,
                           const double beta) const = 0;

   /*!
    * Using a user-specified gradient detection scheme, determine cells which
    * have high gradients and, consequently, should be refined.
    */
   virtual void tagGradientDetectorCells(
                               hier::Patch<DIM>& patch,
                               const double regrid_time,
                               const bool initial_error,
                               const int tag_index,
                               const bool uses_richardson_extrapolation_too)
   {
      NULL_USE(patch);
      NULL_USE(regrid_time);
      NULL_USE(initial_error);
      NULL_USE(tag_index);
      NULL_USE(uses_richardson_extrapolation_too);
   }

   /*!
    * Set user-defined boundary conditions at the physical domain boundary.
    */
   virtual void setPhysicalBoundaryConditions(hier::Patch<DIM>& patch,
                      const double fill_time,
                      const hier::IntVector<DIM>& ghost_width_to_fill) = 0;

   /*!
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual hier::IntVector<DIM> getRefineOpStencilWidth() const
   {
      return(hier::IntVector<DIM>(0));
   }

   /*!
    * Pre- and post-processing routines for implementing user-defined
    * spatial interpolation routines applied to variables.  The 
    * interpolation routines are used in the MOL AMR algorithm
    * for filling patch ghost cells before advancing data on a level
    * and after regridding a level to fill portions of the new level
    * from some coarser level.  These routines are called automatically
    * from within patch boundary filling schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch model.  However, the routines only need to perform some 
    * operations when "USER_DEFINED_REFINE" is given as the interpolation 
    * method for some variable when the patch model registers variables
    * with the MOL integration algorithm, typically.  If the 
    * user does not provide operations that refine such variables in either 
    * of these routines, then they will not be refined.
    *
    * The order in which these operations are used in each patch 
    * boundary filling schedule is:
    * 
    * - \b (1) {Call user's preprocessRefine() routine.}
    * - \b (2) {Refine all variables with standard interpolation operators.}
    * - \b (3) {Call user's postprocessRefine() routine.}
    * 
    * 
    * Also, user routines that implement these functions must use 
    * data corresponding to the d_scratch context on both coarse and
    * fine patches.
    */
   virtual void preprocessRefine(hier::Patch<DIM>& fine,
                                 const hier::Patch<DIM>& coarse,
                                 const hier::Box<DIM>& fine_box,
                                 const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   ///
   virtual void postprocessRefine(hier::Patch<DIM>& fine,
                                  const hier::Patch<DIM>& coarse,
                                  const hier::Box<DIM>& fine_box,
                                  const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /*!
    * Return maximum stencil width needed for user-defined
    * data coarsen operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual hier::IntVector<DIM> getCoarsenOpStencilWidth() const
   {
      return(hier::IntVector<DIM>(0));
   } 

   /*!
    * Pre- and post-processing routines for implementing user-defined
    * spatial coarsening routines applied to variables.  The coarsening 
    * routines are used in the MOL AMR algorithm synchronizing 
    * coarse and fine levels when they have been integrated to the same
    * point.  These routines are called automatically from within the 
    * data synchronization coarsen schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch model.  However, the routines only need to perform some
    * operations when "USER_DEFINED_COARSEN" is given as the coarsening
    * method for some variable when the patch model registers variables
    * with the MOL level integration algorithm, typically.  If the
    * user does not provide operations that coarsen such variables in either
    * of these routines, then they will not be coarsened.
    *
    * The order in which these operations are used in each coarsening
    * schedule is:
    * 
    * - \b (1) {Call user's preprocessCoarsen() routine.}
    * - \b (2) {Coarsen all variables with standard coarsening operators.}
    * - \b (3) {Call user's postprocessCoarsen() routine.}
    * 
    *
    * Also, user routines that implement these functions must use
    * corresponding to the d_new context on both coarse and fine patches
    * for time-dependent quantities.
    */
   virtual void preprocessCoarsen(hier::Patch<DIM>& coarse,
                                  const hier::Patch<DIM>& fine,
                                  const hier::Box<DIM>& coarse_box,
                                  const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   ///
   virtual void postprocessCoarsen(hier::Patch<DIM>& coarse,
                                   const hier::Patch<DIM>& fine,
                                   const hier::Box<DIM>& coarse_box,
                                   const hier::IntVector<DIM>& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   /*!
    * The method of lines integrator controls the context for the data to 
    * be used in the numerical routines implemented in the concrete patch 
    * strategy. These accessor methods allow the patch strategy to access
    * the particular data contexts used in the integrator.
    *
    * Return pointer to data context with ghost cells.
    */
   tbox::Pointer<hier::VariableContext> getInteriorWithGhostsContext() const
   {
      return(d_interior_with_ghosts);
   }

   /*!
    * Return pointer to data context with NO ghosts.
    */
   tbox::Pointer<hier::VariableContext> getInteriorContext() const
   {
      return(d_interior);
   }

   /*!
    * Set pointer to data context with ghosts.
    */
   void setInteriorWithGhostsContext(
      tbox::Pointer<hier::VariableContext> context)
   {
      d_interior_with_ghosts = context;
   }

   /*!
    * Set pointer to data context with NO ghosts.
    */
   void setInteriorContext(
      tbox::Pointer<hier::VariableContext> context)
   {
      d_interior = context;
   }

private:

   tbox::Pointer<hier::VariableContext> d_interior_with_ghosts; 
   tbox::Pointer<hier::VariableContext> d_interior; 
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MethodOfLinesPatchStrategy.C"
#endif
