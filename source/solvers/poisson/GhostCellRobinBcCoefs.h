/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/poisson/GhostCellRobinBcCoefs.h $
 * Package:     SAMRAI solver package
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1818 $
 * Modified:    $LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
 * Description: Robin boundary condition problem-dependent interfaces
 */

#ifndef included_solv_GhostCellRobinBcCoefs
#define included_solv_GhostCellRobinBcCoefs


#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif


/*
 * SAMRAI classes
 */

#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif

#ifndef included_hier_Patch
#include "Patch.h"
#endif

#ifndef included_pdat_ArrayData
#include "ArrayData.h"
#endif

#ifndef included_solv_RobinBcCoefStrategy
#include "RobinBcCoefStrategy.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif


namespace SAMRAI {
    namespace solv {


/*!
 * @brief A prefabricated Robin boundary condition coefficients
 * for the case where cell-centered data is fixed at the first
 * ghost-cell centers.
 *
 * This class is intended to make the use of the Robin boundary
 * condition utterly trivial for users who who already have the
 * correct values set in the ghost cells.  The motivation for
 * this seemingly needless task is to interpret the requirement
 * that solution is fixed at ghost cell centers to solvers that
 * do not operate directly on the ghost cell values.
 * An example is linear solvers that operate on Ax=b and
 * require that the boundary condition be written as changes
 * to A and b.
 *
 * This implementation of the strategy
 * class RobinBcCoefStrategy<DIM> can be used when ghost cell
 * values are known and have been written to the ghost cells
 * of the data being set.  You provide the patch data index
 * to the cell-centered data, defined with a non-zero ghost
 * cell width, where the ghost cell values can be found.
 *
 * This implementation corresponds to a specific discretization
 * of the Robin formula described in RobinBcCoefStrategy<DIM>.
 * It assumes a linear variation of the data between the first
 * interior and first ghost cells.  It sets up the coefficients
 * such that this linear extrapolation gives the correct value
 * at the ghost cell center.  This results in the coefficient
 * a being 1.0/(1+0.5*h) and g being the a times the ghost cell
 * value.  h is the grid spacing normal to the boundary.
 */
template<int DIM> class GhostCellRobinBcCoefs
  : public RobinBcCoefStrategy<DIM>
{

public:

   /*!
    * @brief Constructor
    *
    * @param object_name Name of object for output purposes.
    */
   GhostCellRobinBcCoefs( std::string object_name="" );


   /*!
    * @brief Destructor.
    */
   virtual ~GhostCellRobinBcCoefs<DIM>(void);


   /*!
    * @brief Function to fill arrays of Robin boundary
    * condition coefficients at a patch boundary.
    *
    * This implementation of the virtual function
    * RobinBcCoefStrategy<DIM>::setBcCoefs()
    * sets up the coefficients as described in
    * the above notes.
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
    * @param fill_time Solution time corresponding to filling,
    *        for use when coefficients are time-dependent.
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


   /*!
    * @brief Set the patch data index of the data providing
    * ghost cell values.
    *
    * The index must correspond to cell-centered double
    * data with the given ghost width.
    *
    * @param ghost_data_id patch data index of ghost data
    * @param extensions_fillable the number of extensions past
    *        edge of a patch that has valid ghost cell values.
    */
   void setGhostDataId( int ghost_data_id,
                        hier::IntVector<DIM> extensions_fillable= hier::IntVector<DIM>(0) );


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
    * @brief Object name.
    */
   std::string d_object_name;

   /*
    * @brief hier::Index of cell-centered double
    * data to provide ghost cell values.
    *
    * Set to -1 until setGhostDataId() is called.
    */
   int d_ghost_data_id;

   /*
    * @brief Extensions fillable that can be used for data at d_ghost_data_id.
    */
   hier::IntVector<DIM> d_extensions_fillable;

   /*
    * @brief tbox::Timer classes for performance measurement
    */
   tbox::Pointer<tbox::Timer> t_set_bc_coefs;

};

}
}

#endif	// included_solv_GhostCellRobinBcCoefs

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "GhostCellRobinBcCoefs.C"
#endif
