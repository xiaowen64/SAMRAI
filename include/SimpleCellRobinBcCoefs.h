/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/SimpleCellRobinBcCoefs.h $
 * Package:     SAMRAI solver package
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2132 $
 * Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
 * Description: Robin boundary condition problem-dependent interfaces
 */

#ifndef included_solv_SimpleCellRobinBcCoefs
#define included_solv_SimpleCellRobinBcCoefs


#include "SAMRAI_config.h"


/*
 * SAMRAI classes
 */

#include "BoundaryBox.h"

#include "PatchHierarchy.h"

#include "Patch.h"

#include "ArrayData.h"

#include "RobinBcCoefStrategy.h"

#include "tbox/Pointer.h"


namespace SAMRAI {
    namespace solv {


/*!
 * @brief A prefabricated Robin boundary condition coefficients
 * with an interface like the old Poisson solvers.
 *
 * This class is intended to make it easy for users of the old
 * Poisson solver to adapt to the new solver by providing an
 * interface similar to that of the old solver.  Underneath,
 * the boundary condition is converted to Robin bc coefficients
 * used by the new solver.
 *
 * This class refers to some grid-based outerside or ghost cell
 * data touching patch outer sides when providing the Robin bc
 * coefficients.  In the most general case, it is unable to
 * provide coefficients for patches outside the grid on which
 * that data is provided.  It is also unable to provide coefficients
 * for boundary boxes that extend past the edge or corner of a
 * patch.  This may limit this class from being used when certain
 * features of the Poisson solver is enabled.
 */
template<int DIM> class SimpleCellRobinBcCoefs
  : public RobinBcCoefStrategy<DIM>
{

public:

   /*!
    * @brief Constructor
    */
   SimpleCellRobinBcCoefs(const std::string& object_name=std::string());


   /*!
    * @brief Destructor.
    */
   virtual ~SimpleCellRobinBcCoefs<DIM>(void);



   //@{ @name Inherited from RobinBcCoefStrategy<DIM>

   /*!
    * @brief Function to fill arrays of Robin boundary
    * condition coefficients at a patch boundary.
    *
    * This implementation of the virtual function
    * RobinBcCoefStrategy<DIM>::setBcCoefs()
    * uses information specified by the call to setBoundaries()
    * to determine the coefficients.
    *
    * @param acoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    *        If this is a null pointer, then the calling function
    *        is not interested in a, and you can disregard it.
    * @param bcoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    * @param gcoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    * @param variable variable to set the coefficients for.
    * @param patch patch requiring bc coefficients
    * @param bdry_box boundary box showing where on the boundary
    *        the coefficient data is needed.
    * @param fill_time solution time corresponding to filling, for use
    *        when coefficients are time-dependent.
    */
   void setBcCoefs (
      tbox::Pointer<pdat::ArrayData<DIM,double> > &acoef_data ,
      tbox::Pointer<pdat::ArrayData<DIM,double> > &bcoef_data ,
      tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
      const tbox::Pointer< hier::Variable<DIM> > &variable ,
      const hier::Patch<DIM> &patch ,
      const hier::BoundaryBox<DIM> &bdry_box ,
      double fill_time=0.0 ) const;

   hier::IntVector<DIM> numberOfExtensionsFillable() const;

   //@}




   /*!
    * @brief Set the hierarchy where boundary data associated
    * with the hierarchy is found.
    *
    * This class requires you to specify some grid data
    * associated with a hierarchy, such as the Dirichlet
    * boundary values, the flux or the Dirichlet/Neumann
    * flag.  That hierarchy and the range of relevant
    * patch levels is specified by calling this function.
    */
   void setHierarchy(tbox::Pointer< hier::PatchHierarchy<DIM> >,
                     const int ln_min=-1,
                     const int ln_max=-1);


   /*!
    * @brief Specify the boundary conditions that are to be used at the
    * physical domain boundary.
    *
    * The boundary conditions specified as the
    * string argument "boundary_type."  The boundary type argument can be
    * "Dirichlet", "Neumann", or "Mixed".
    *
    * If using Dirichlet boundary conditions, then before the solver is
    * called, the storage for the unknown u
    * must have a layer of ghost cells at least one cell wide that includes
    * the Dirichlet boundary values.
    *
    * If using Neumann boundary conditions, then before the solver is called,
    * the outerface boundary flux data must be set for the Neumann conditions.
    * The fluxes argument gives the patch data index of this flux
    * data.
    *
    * The mixed boundary type is for a mixture of Dirichlet and Neumann
    * boundary conditions are used at the physical domain boundary.
    * The fluxes argument gives the patch data index of the outerface data
    * that specifies the flux data for the Neumann conditions.  The flags
    * array is an outerface data array of integer flags that specifies whether
    * Dirichlet (flag == zero) or Neumann (flag == one) conditions are to be
    * used at a particular cell boundary face.  Note that the flag data must
    * be set before the matrix entries can be computed and the flux data
    * must be set before the solver is called.  The bdry_types argument can
    * be used if the boundary conditions are mixed but one or more of the
    * faces of the physical boundary are entirely either Dirichlet or
    * Neumann boundaries.  The bdry_types argument should be an array of
    * 2*DIM integers, specifying the boundary conditions on each side of
    * the physical domain.  It should be ordered {x_lo, x_hi, y_lo, y_hi,
    * z_lo, z_hi}, with the values for each face being 0 for Dirichlet
    * conditions, 1 for Neumann conditions, and 2 for mixed boundary
    * conditions.  The bdry_type argument is never required, but if used
    * it can sometimes make the PoissonHYPRESolver class more efficient.
    */

   void setBoundaries(const std::string& boundary_type,
                      const int fluxes = -1,
                      const int flags = -1,
                      int* bdry_types = NULL);


   /*!
    * @brief Cache data providing Dirichlet boundary values.
    *
    * This function makes a private copy of the relevant ghost cell data
    * that it later uses provide the coefficient g on Dirichlet boundaries.
    * The index must correspond to cell-centered double
    * data with non-zero ghost width.
    *
    * Functions setHierarchy() and setBoundaries()
    * should be called before this one.
    * This function should be called each time the hierarchy or
    * Dirichlet data changes.
    *
    * @param dirichlet_data_id patch data id of the source cell data
    *        for copy.
    */
   void cacheDirichletData( int dirichlet_data_id );


   /*!
    * @brief Copy cached Dirichlet data into ghost cells.
    *
    * Reverse action of cacheDirichletData by copying cached data back
    * into the ghost cells.
    *
    * The cached data is not dallocated.
    *
    * @param dirichlet_data_id patch data id of the destination cell data
    *        for copy.
    */
   void restoreDirichletData( int dirichlet_data_id );


   /*!
    * @brief Set the patch data index of the diffusion
    * coefficient used in Neumann boundary fluxes.
    *
    * The diffusion coefficient, along with the prescribed flux,
    * is used to set the gradient of the solution normal to the boundary.
    * By default, the diffusion coefficient is assumed to be 1.
    * If used, the diffusion coefficient data id must be set before
    * asking for the coefficient g, which depends on it.
    * The index must correspond to side-centered double data.
    *
    * This function overrides the effect of setDiffusionCoefConstant().
    */
   void setDiffusionCoefId( int diffusion_coef_id );


   /*!
    * @brief Set the value of the diffusion coefficient
    * used in Neumann boundary fluxes to a constant.
    *
    * This function is similar to setDiffusionCoefId() but is
    * used when the diffusion coefficient is a constant.
    *
    * This function overrides the effect of setDiffusionCoefId().
    */
   void setDiffusionCoefConstant( double diffusion_coef_value );


private:

   /*!
    * @brief Return box describing the index space of surfaces
    * defined by a boundary box.
    *
    * Define a box describing the indices of the surfaces corresponding
    * to the input boundary box.  A surface is a face in 3D and an edge
    * in 2D.  These surfaces lie on the boundary itself.
    *
    * The input boundary_box must be of type 1
    * (see hier::BoundaryBox::getBoundaryType()).
    *
    * This is a utility function for working with the surface
    * indices coresponding to a boundary box.
    *
    * @param boundary_box input boundary box
    * @return a box to define the face indices corresponding to
    *    boundary_box
    */
   hier::Box<DIM> makeSideBoundaryBox(
      const hier::BoundaryBox<DIM> &boundary_box ) const;

   /*!
    * @brief object name
    */
   std::string d_object_name;
   tbox::Pointer< hier::PatchHierarchy<DIM> > d_hierarchy;
   int d_ln_min;
   int d_ln_max;
   /*!
    * @brief array of boundary type on each side
    */
   int d_bdry_types[2*DIM];
   /*!
    * @brief patch index for fluxes
    */
   int d_flux_id;
   /*!
    * @brief patch index for flags
    */
   int d_flag_id;
   /*!
    * @brief patch index for Dirichlet values.
    */
   int d_dirichlet_data_id;
   /*!
    * @brief patch index for diffusion coefficients if it is variable.
    */
   int d_diffusion_coef_id;
   /*!
    * @brief value of for diffusion coefficients if it is constant.
    */
   double d_diffusion_coef_constant;
   /*!
    * @brief Cached ghost cell value used in Dirichlet bc.
    *
    * Cached boundary box ghost cell data are stored in this 1D
    * array.  For the position of a particular box, see
    * d_dirichlet_data_position.
    */
   tbox::Array<tbox::Pointer<pdat::ArrayData<DIM,double> > > d_dirichlet_data;
   /*!
    * @brief Position of cached boundary boxes of ghost cell data.
    *
    * The position of the cached boundary box bn of patch pn of
    * level ln is d_dirichlet_data_pos[ln][pn]+bn.
    */
   tbox::Array<tbox::Array<int> > d_dirichlet_data_pos;

   /*!
    * @brief Timers for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_set_bc_coefs;

};

}
}

#endif	// included_solv_SimpleCellRobinBcCoefs

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SimpleCellRobinBcCoefs.C"
#endif
