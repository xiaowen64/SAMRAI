#ifndef included_StencilXD
#define included_StencilXD

#include <string>

#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"

using namespace SAMRAI;

class Stencil :
  public algs::HyperbolicPatchStrategy,
  public appu::BoundaryUtilityStrategy
{
  public:
    Stencil(
        const std::string& name,
        const tbox::Dimension& dim,
        boost::shared_ptr<tbox::Database> input_db,
        boost::shared_ptr<geom::CartesianGridGeometry> grid_geom);

    void
      registerModelVariables(
          algs::HyperbolicLevelIntegrator* integrator);

    void
      initializeDataOnPatch(
          hier::Patch& patch,
          const double data_time,
          const bool initial_time);

    double
      computeStableDtOnPatch(
          hier::Patch& patch,
          const bool initial_time,
          const double dt_time);

    void
      computeFluxesOnPatch(
          hier::Patch& patch,
          const double time,
          const double dt);

    void
      conservativeDifferenceOnPatch(
          hier::Patch& patch,
          const double time,
          const double dt,
          bool at_syncronization);

    void
      tagGradientDetectorCells(
          hier::Patch& patch,
          const double regrid_time,
          const bool initial_error,
          const int tag_index,
          const bool uses_richardson_extrapolation_too);

    void
      setPhysicalBoundaryConditions(
          hier::Patch& patch,
          const double fill_time,
          const hier::IntVector& ghost_width_to_fill);

   /**
    * Return stencil width of conservative linear interpolation operations.
    */
   hier::IntVector
   getRefineOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector::getOne(dim);
   }

   void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /**
    * Refine velocity and pressure from coarse patch to fine patch
    * so that momentum and total energy are conserved.
    */
   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   /**
    * Return stencil width of conservative averaging operations.
    */
   hier::IntVector
   getCoarsenOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector::getZero(dim);
   }

   void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   /**
    * Coarsen velocity and pressure from coarse patch to fine patch
    * so that momentum and total energy are conserved.
    */
   void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   void
   readDirichletBoundaryDataEntry(
      const boost::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It is a blank implementation
    * for the purposes of this class.
    */
   void
   readNeumannBoundaryDataEntry(
      const boost::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   void
     registerVisItDataWriter(
        boost::shared_ptr<appu::VisItDataWriter> viz_writer);

  private:
    std::string d_object_name;

    boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

    const tbox::Dimension d_dim;

    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_rho_variables;

    hier::IntVector d_nghosts;

    boost::shared_ptr<appu::VisItDataWriter> d_visit_writer;
};

#endif
