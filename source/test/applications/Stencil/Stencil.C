#include "Stencil.h"

#include "SAMRAI/pdat/CellData.h"

#include "SAMRAI/tbox/RAJA_API.h"

template <typename T, int DIM>
struct CellView : 
  public tbox::ArrayView<DIM, T>
{
  SAMRAI_INLINE CellView(const boost::shared_ptr<pdat::CellData<T> >& data, int depth = 0)
    : tbox::ArrayView<DIM, T>(data->getPointer(depth), data->getGhostBox(), depth)
  {
  }
};


Stencil::Stencil(
        const std::string& name,
        const tbox::Dimension& dim,
        boost::shared_ptr<tbox::Database> input_db,
        boost::shared_ptr<geom::CartesianGridGeometry> grid_geom):
  algs::HyperbolicPatchStrategy(),
  d_object_name(name),
  d_grid_geometry(grid_geom),
  d_dim(dim),
  d_rho_variables(),
  d_nghosts(hier::IntVector(dim, 2))
{
  const int num_variables = input_db.getIntegerWithDefault("num_variables", 1);

  for (int i = 0; i < num_variables; ++i) {
    std::ostringstream oss;
    oss << "rho_" << i;
    std::string var_name = oss.str();

    d_rho_variables.push_back(new pdat::CellVariable<double>(dim, var_name, 1));
  }
}

void
Stencil::registerModelVariables(
  algs::HyperbolicLevelIntegrator* integrator)
{
  for ( const boost::shared_ptr<pdat::CellVariable<double> > rho_var : d_rho_variables ) {
    integrator->registerVariable(
        rho_var, 
        d_nghosts,
        algs::HyperbolicLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
  }
}

void
Stencil::initializeDataOnPatch(
  hier::Patch& patch,
  const double data_time,
  const bool initial_time)
{
  // initialize
  for ( const boost::shared_ptr<pdat::CellVariable<double> > rho_var : d_rho_variables ) {
    CellView<double, 2> rho(BOOST_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext())));

    tbox::for_all2<tbox::policy::parallel>(patch.getBox(), [=] __host__ __device__ (int k, int j) {
      rho(j,k) = 0.0;
    });
  }
}

double
Stencil::computeStableDtOnPatch(
  hier::Patch& patch,
  const bool initial_time,
  const double dt_time)
{

  // TODO: calc dt

  return 0.5;
}

void
Stencil::computeFluxesOnPatch(
  hier::Patch& patch,
  const double time,
  const double dt)
{
  // TODO: do some work
}

void
Stencil::conservativeDifferenceOnPatch(
  hier::Patch& patch,
  const double time,
  const double dt,
  bool at_syncronization)
{
  // This function is always a no-op
  (void) patch;
}

void
Stencil::tagGradientDetectorCells(
  hier::Patch& patch,
  const double regrid_time,
  const bool initial_error,
  const int tag_index,
  const bool uses_richardson_extrapolation_too)
{
  /*
   * Only need to tag the first variable.
   */
  CellView<double, 2> rho(BOOST_CAST<pdat::CellData<double> >(patch.getPatchData(d_rho_variables[0], getDataContext())));

  CellView<int, 2> tags(BOOST_CAST<pdat::CellData<int> >(patch.getPatchData(tag_index)));

  Real tag_threshold = d_tag_threshold;

  tbox::for_all2<tbox::policy::parallel>(patch.getBox(), [=] __device__ (int k, int j) {
    Real d2x = ABS(rho(j+1,k) - 2.0*rho(j,k) + rho(j-1,k));
    Real d2y = ABS(rho(j,k+1) - 2.0*rho(j,k) + rho(j,k-1));

    Real dxy = ABS(rho(j+1,k+1) - 2.0*rho(j,k) + rho(j-1,k-1));
    Real dyx = ABS(rho(j-1,k+1) - 2.0*rho(j,k) + rho(j+1,k-1));

    Real dd = MAX(d2x,MAX(d2y,MAX(dxy,dyx)));

    if (dd > tag_threshold) {
      tags(j,k) = 1;
    }
  });
}

void
Stencil::setPhysicalBoundaryConditions(
  hier::Patch& patch,
  const double fill_time,
  const hier::IntVector& ghost_width_to_fill)
{
  const int depth = ghost_width_to_fill[0];

  const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SHARED_PTR_CAST(geom::CartesianPatchGeometry,
        patch.getPatchGeometry()));

  const std::vector<hier::BoundaryBox>& edge_bdry 
    = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);


  for ( const boost::shared_ptr<pdat::CellVariable<double> > rho_var : d_rho_variables ) {
    CellView<double, 2> field(BOOST_CAST<pdat::CellData<double> >(patch.getPatchData(rho_var, getDataContext())));

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    const int ifirst0 = ifirst(0);
    const int ifirst1 = ifirst(1);
    const int ilast0 = ilast(0);
    const int ilast1 = ilast(1);

    for(int i = 0; i < edge_bdry.size(); i++) {

      auto edge = edge_bdry[i].getLocationIndex();

      hier::Box boundary_box(pgeom->getBoundaryFillBox(
            edge_bdry[i],
            patch.getBox(),
            ghost_width_to_fill);

          switch(edge) {
          case (BdryLoc::YLO) :
          {
          tbox::for_all1<tbox::policy::parallel>(
              boundary_box,
              0,
              [=] __host__ __device__ (int j) {
              field(j, ifirst1-1)  = field(j,ifirst1);
              if (depth == 2) {
              field(j, ifirst1-2)  = field(j,ifirst1+1);
              }
              });
          }
          break;
          case (BdryLoc::YHI) :
          {
            tbox::for_all1<tbox::policy::parallel>(
                boundary_box,
                0,
                [=] __host__ __device__ (int j) {
                field(j,ilast1+1) = field(j,ilast1);
                if (depth == 2) {
                field(j,ilast1+2) = field(j,ilast1-1);
                }
                });
          }
          break;
          case (BdryLoc::XLO) :
          {
            tbox::for_all1<tbox::policy::parallel>(
                boundary_box,
                1,
                [=] __host__ __device__ (int k) {
                field(ifirst0-1,k) = field(ifirst0,k);
                if (depth == 2) {
                field(ifirst0-2,k) = field(ifirst0+1,k);
                }
                });
          }
          break;
          case (BdryLoc::XHI) :
          {
            tbox::for_all1<tbox::policy::parallel>(
                boundary_box,
                1,
                [=] __host__ __device__ (int k) {
                field(ilast0+1,k) = field(ilast0,k);
                if (depth == 2) {
                field(ilast0+2,k) = field(ilast0-1,k);
                }
                });
          }
          break;
          }
    }
  }
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
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
  // no-op
}

 void
 Stencil::readNeumannBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
  // no-op
}
