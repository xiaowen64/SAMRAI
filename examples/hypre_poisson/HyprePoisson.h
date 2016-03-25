/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/hypre_poisson/HyprePoisson.h $
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Example user class for solving Poisson using Hypre.
 */

#ifndef included_HyprePoisson
#define included_HyprePoisson

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

using namespace std;

#if !defined(HAVE_HYPRE)

/*
*************************************************************************
* If the library is not compiled with hypre, print an error.
* If we're running autotests, skip the error and compile an empty
* class.
*************************************************************************
*/
#if (TESTING != 1)
#error "This example requires SAMRAI be compiled with hypre."
#endif

#else

#include "CartesianVizamraiDataWriter.h"
#include "CellVariable.h"
#include "tbox/Database.h"
#include "Box.h"
#include "LocationIndexRobinBcCoefs.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "IntVector.h"
#include "CellPoissonHypreSolver.h"
#include "SideVariable.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "VisDerivedDataStrategy.h"
#include "VisItDataWriter.h"

namespace SAMRAI {

/*!
 * @brief Class to solve a sample Poisson equation on a SAMR grid.
 *
 * This class demonstrates how use the HYPRE Poisson solver
 * class to solve Poisson's equation on a single level
 * within a hierarchy.
 *
 * We set up and solve the following problem:
 *
 *   2d: div(grad(u)) = -2 (pi^2) sin(pi x) sin(pi y)
 *
 *   3d: div(grad(u)) = -3 (pi^2) sin(pi x) sin(pi y) sin(pi z)
 *
 * which has the exact solution
 *
 *   2d: u = sin(pi x) sin(pi y)
 *
 *   3d: u = sin(pi x) sin(pi y) sin(pi z)
 *
 * This class inherits and implements virtual functions from
 * - mesh::StandardTagAndInitStrategy<NDIM> to initialize data
 *   on the SAMR grid.
 * - appu::VisDerivedDataStrategy<NDIM> to write out certain data
 *   in a vis file, such as the error of the solution.
 */
class HyprePoisson :
   public mesh::StandardTagAndInitStrategy<NDIM> ,
   public appu::VisDerivedDataStrategy<NDIM>
{

public:

   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    *
    * @param object_name Ojbect name
    */
   HyprePoisson( const string &object_name,
                 tbox::Pointer<tbox::Database> database=NULL);

   virtual ~HyprePoisson();



   //@{ @name mesh::StandardTagAndInitStrategy<NDIM> virtuals


   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on
    * the grid.  All the information needed to initialize the grid
    * are in the arguments.
    *
    * @see mesh::StandardTagAndInitStrategy<NDIM>::initializeLevelData()
    */
   virtual void initializeLevelData (
      const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy ,
      const int level_number ,
      const double init_data_time ,
      const bool can_be_refined ,
      const bool initial_time ,
      const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level = 
                       tbox::Pointer<hier::BasePatchLevel<NDIM> >((0)) ,
      const bool allocate_data = true );

   /*!
    * @brief Reset any internal hierarchy-dependent information.
    */
   virtual void resetHierarchyConfiguration (
      tbox::Pointer<hier::BasePatchHierarchy<NDIM> > new_hierarchy ,
      int coarsest_level ,
      int finest_level );

   //@}


   //@{ @name appu::VisDerivedDataStrategy<NDIM> virtuals

   virtual bool packDerivedDataIntoDoubleBuffer(
      double *buffer ,
      const hier::Patch<NDIM> &patch ,
      const hier::Box<NDIM> &region ,
      const string &variable_name ,
      int depth_id ) const;

   //@}



   /*!
    * @brief Solve using HYPRE Poisson solver
    *
    * Set up the linear algebra problem and use a
    * solv::CellPoissonHypreSolver<NDIM> object to solve it.
    * -# Set initial guess
    * -# Set boundary conditions
    * -# Specify Poisson equation parameters
    * -# Call solver
    *
    * @return whether solver converged
    */
   bool solvePoisson();

   /*!
    * @brief Set up external plotter to plot internal
    * data from this class.
    *
    * After calling this function, the external
    * data writer may be used to write the
    * viz file for this object.
    *
    * The internal hierarchy is used and must be
    * established before calling this function.
    * (This is commonly done by building a hierarchy
    * with the mesh::StandardTagAndInitStrategy<NDIM> virtual
    * functions implemented by this class.)
    *
    * @param viz_writer Vizramrai or VisIt data writer
    */
   int setupExternalPlotter(
      appu::CartesianVizamraiDataWriter<NDIM> &viz_writer ) const;

#ifdef HAVE_HDF5
   int setupExternalPlotter(
      appu::VisItDataWriter<NDIM> &viz_writer ) const;
#endif

private:

  string d_object_name;

  tbox::Pointer<hier::PatchHierarchy<NDIM> > d_hierarchy;

  //@{
  /*!
   * @name Major algorithm objects.
   */

  /*!
   * @brief HYPRE poisson solver.
   */
  solv::CellPoissonHypreSolver<NDIM> d_poisson_hypre;

  /*!
   * @brief Boundary condition coefficient implementation.
   */
  solv::LocationIndexRobinBcCoefs<NDIM> d_bc_coefs;

  //@}

  //@{
private:
  /*!
   * @name Private state variables for solution.
   */

  /*!
   * @brief Context owned by this object.
   */
  tbox::Pointer<hier::VariableContext> d_context;

  /*!
   * @brief Descriptor indices of internal data.
   *
   * These are initialized in the constructor and never change.
   */
  int d_comp_soln_id, d_exact_id, d_rhs_id;

  //@}

};


}

#endif
#endif	// included_HyprePoisson
