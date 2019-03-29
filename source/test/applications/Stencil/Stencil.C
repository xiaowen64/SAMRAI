#include "Stencil.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/tbox/Collectives.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "SAMRAI/tbox/NVTXUtilities.h"

#include "math.h"

#define MAX(a, b) fmax(a, b)
#define MIN(a, b) fmin(a, b)
#define ABS(a) fabs(a)

Stencil::Stencil(
   const std::string& name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<geom::CartesianGridGeometry> grid_geom):
   algs::HyperbolicPatchStrategy(),
   d_object_name(name),
   d_grid_geometry(grid_geom),
   d_velocity({1.0,0.0}),
   d_dim(dim),
   d_rho_variables(),
   d_rho_update(),
   d_nghosts(hier::IntVector(dim, 2))
{
   const int num_variables = input_db->getIntegerWithDefault("num_variables", 1);

   d_tag_threshold = input_db->getDoubleWithDefault("tag_threshold", 0.5);

   d_rho_update = std::make_shared<pdat::CellVariable<double> >(dim, "update", 1);

   for (int i = 0; i < num_variables; ++i) {
      std::ostringstream oss;
      oss << "rho_" << i;
      std::string var_name = oss.str();

      d_rho_variables.push_back( std::make_shared<pdat::CellVariable<double> >(dim, var_name, 1));
   }
}

void
Stencil::registerModelVariables(
   algs::HyperbolicLevelIntegrator* integrator)
{
   integrator->registerVariable(
      d_rho_update,
      hier::IntVector::getZero(d_dim),
      algs::HyperbolicLevelIntegrator::TEMPORARY,
      d_grid_geometry);

   int i = 0;
   for ( const auto& rho_var : d_rho_variables ) {

      integrator->registerVariable(
         rho_var,
         d_nghosts,
         algs::HyperbolicLevelIntegrator::TIME_DEP,
         d_grid_geometry,
         "CONSERVATIVE_COARSEN",
         "CONSERVATIVE_LINEAR_REFINE");

      if (i == 0) {
         hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

         d_visit_writer->registerPlotQuantity(
            rho_var->getName(),
            "SCALAR",
            vardb->mapVariableAndContextToIndex(
               rho_var, integrator->getPlotContext()));
      }
      i++;
   }
}

void
Stencil::initializeDataOnPatch(
   hier::Patch& patch,
   const double data_time,
   const bool initial_time)
{
   RANGE_PUSH("Stencil::init", 1);

   // initialize
   if (initial_time) {
      for ( const auto& rho_var : d_rho_variables ) {
         auto rho = pdat::get_view<2, pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext()));
         // CellView<double, 2> rho(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext())));

         pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
            rho(j,k) = 0.0;
         });
      }
   }

   RANGE_POP;
}

double
Stencil::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time,
   const double dt_time)
{
   RANGE_PUSH("Stencil::dt", 1);

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(patch.getPatchGeometry()));

   const double dx = pgeom->getDx()[0];
   const double dy = pgeom->getDx()[1];
   const double min_length = dx > dy ? dy : dx;

   double velocity_mag = std::sqrt(d_velocity[0]*d_velocity[0] + d_velocity[1]*d_velocity[1]);

   RANGE_POP;

   return min_length / velocity_mag;
}

void
Stencil::computeFluxesOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt)
{ }

void
Stencil::conservativeDifferenceOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt,
   bool at_syncronization)
{
   RANGE_PUSH("Stencil::conservativeDifference", 1);

   auto rho_new = pdat::get_view<2, pdat::CellData<double>>(patch.getPatchData(d_rho_update, getDataContext()));
   // CellView<double, 2> rhoNew(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(d_rho_update, getDataContext())));

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(patch.getPatchGeometry()));

   const double dx = pgeom->getDx()[0];
   const double dy = pgeom->getDx()[1];

   for ( const auto& rho_var : d_rho_variables ) {
      auto rho = pdat::get_view<2, pdat::CellData<double>>(patch.getPatchData(rho_var, getDataContext()));
      // CellView<double, 2> rho(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext())));

      const double an_abs_lr = ABS(d_velocity[0]) * 0.5;
      const double an_abs_tb = ABS(d_velocity[1]) * 0.5;
      const double an_lr = d_velocity[0] * 0.5;
      const double an_tb = d_velocity[1] * 0.5;

      pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
         const double FL = (rho(j-1,k) + rho(j,k)) * an_abs_lr + (rho(j-1,k) - rho(j,k)) * an_lr;
         const double FR = (rho(j,k) + rho(j+1,k)) * an_abs_lr + (rho(j,k) - rho(j+1,k)) * an_lr;
         const double FT = (rho(j,k) + rho(j,k+1)) * an_abs_tb + (rho(j,k) - rho(j,k+1)) * an_tb;
         const double FB = (rho(j,k-1) + rho(j,k)) * an_abs_tb + (rho(j,k-1) - rho(j,k)) * an_tb;
         rho_new(j,k) = rho(j,k) - dt * (1.0 / dx * (FR - FL) + 1.0 / dy * (FT - FB));
      });

      pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
         rho(j,k) = rho_new(j,k);
      });
   }

   RANGE_POP;
}

void
Stencil::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   RANGE_PUSH("Stencil::tag", 1);
   /*
    * Only need to tag the first variable.
    */
   auto rho = pdat::get_const_view<2, pdat::CellData<double> >(patch.getPatchData(d_rho_variables[0], getDataContext()));
   // CellView<double, 2> rho(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(d_rho_variables[0], getDataContext())));

   auto tags = pdat::get_view<2, pdat::CellData<double> >(patch.getPatchData(tag_index));
   // CellView<int, 2> tags(SAMRAI_SHARED_PTR_CAST<pdat::CellData<int> >(patch.getPatchData(tag_index)));

   double tag_threshold = d_tag_threshold;

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   const int ifirst0 = ifirst(0);
   const int ifirst1 = ifirst(1);
   const int ilast0 = ilast(0);
   const int ilast1 = ilast(1);

   pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
      const double d2x = ABS(rho(j+1,k) - 2.0*rho(j,k) + rho(j-1,k));
      const double d2y = ABS(rho(j,k+1) - 2.0*rho(j,k) + rho(j,k-1));

      /*
       * TODO: fix boundary conditions for diagonal gradient detection
       */
      const double dxy = 0.0; // ABS(rho(j+1,k+1) - 2.0*rho(j,k) + rho(j-1,k-1));
      const double dyx = 0.0; // ABS(rho(j-1,k+1) - 2.0*rho(j,k) + rho(j+1,k-1));

      const double dd = MAX(d2x,MAX(d2y,MAX(dxy,dyx)));

      if (dd > tag_threshold) {
         tags(j,k) = 1;
      }
   });

   RANGE_POP;
}

void
Stencil::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   RANGE_PUSH("Stencil::boundaries", 1);

   const int depth = ghost_width_to_fill[0];

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(patch.getPatchGeometry()));

   const std::vector<hier::BoundaryBox>& edge_bdry
      = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);


   for ( const auto& rho_var : d_rho_variables ) {
      auto field = pdat::get_view<2, pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext()));
      // CellView<double, 2> field(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext())));

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();

      const int ifirst0 = ifirst(0);
      const int ifirst1 = ifirst(1);
      const int ilast0 = ilast(0);
      const int ilast1 = ilast(1);

      for(int i = 0; i < edge_bdry.size(); i++) {

         const auto edge = edge_bdry[i].getLocationIndex();

         const hier::Box boundary_box(pgeom->getBoundaryFillBox(
                                         edge_bdry[i],
                                         patch.getBox(),
                                         ghost_width_to_fill));

         switch(edge) {
         case (BdryLoc::YLO) :
         {
            pdat::parallel_for_all(boundary_box, 0, [=] SAMRAI_HOST_DEVICE (int j) {
               field(j, ifirst1-1)  = field(j,ifirst1);
               if (depth == 2) { field(j, ifirst1-2)  = field(j,ifirst1+1); }
            });
         }
         break;
         case (BdryLoc::YHI) :
         {
            pdat::parallel_for_all(boundary_box, 0, [=] SAMRAI_HOST_DEVICE (int j) {
               field(j,ilast1+1) = field(j,ilast1);
               if (depth == 2) { field(j,ilast1+2) = field(j,ilast1-1); }
            });
         }
         break;
         case (BdryLoc::XLO) :
         {
            pdat::parallel_for_all(boundary_box, 1, [=] SAMRAI_HOST_DEVICE (int k) {
               field(ifirst0-1,k) = 1.0;
               if (depth == 2) { field(ifirst0-2,k) = 1.0; }
            });
         }
         break;
         case (BdryLoc::XHI) :
         {
            pdat::parallel_for_all(boundary_box, 1, [=] SAMRAI_HOST_DEVICE (int k) {
               field(ilast0+1,k) = field(ilast0,k);
               if (depth == 2) { field(ilast0+2,k) = field(ilast0-1,k); }
            });
         }
         break;
         }
      }
   } // for rho_var

   RANGE_POP;
}

double
Stencil::computeNorm(const std::shared_ptr<hier::VariableContext>& context, hier::Patch& patch) const
{
   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(patch.getPatchGeometry()));

   const double dx = pgeom->getDx()[0];
   const double dy = pgeom->getDx()[1];

   tbox::parallel_reduction_variable_t<tbox::Reduction::Sum, double> norm(0.0);
   for ( const auto& rho_var : d_rho_variables ) {
      auto rho = pdat::get_const_view<2, pdat::CellData<double>>(patch.getPatchData(rho_var, context));
      // CellView<double, 2> rho(SAMRAI_SHARED_PTR_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, context)));
      pdat::parallel_for_all(patch.getBox(), [=] SAMRAI_HOST_DEVICE (int k, int j) {
         norm += ABS(rho(j,k)) * dx * dy;
      });
   }
   return norm.get();
}

void
Stencil::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   // no-op
}

void
Stencil::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   // no-op
}

void
Stencil::readDirichletBoundaryDataEntry(
   const std::shared_ptr<tbox::Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   // no-op
}

void
Stencil::readNeumannBoundaryDataEntry(
   const std::shared_ptr<tbox::Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   // no-op
}

#if 1
void
Stencil::registerVisItDataWriter(std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   d_visit_writer = viz_writer;
}
#endif
