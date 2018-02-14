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
  d_rho(new pdat::CellVariable<double>(dim, "rho", 1)),
  d_nghosts(hier::IntVector(dim, 2))
{
}


void
Stencil::registerModelVariables(
  algs::HyperbolicLevelIntegrator* integrator)
{
  integrator->registerVariable(
      d_rho, 
      d_nghosts,
      algs::HyperbolicLevelIntegrator::TIME_DEP,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "CONSERVATIVE_LINEAR_REFINE");

  hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

  d_visit_writer->registerPlotQuantity(
    "rho",
     "SCALAR",
     vardb->mapVariableAndContextToIndex(
        d_rho, integrator->getPlotContext()));
}

void
Stencil::initializeDataOnPatch(
  hier::Patch& patch,
  const double data_time,
  const bool initial_time)
{
  // initialize
  
  CellView<double, 2> rho(BOOST_CAST<pdat::CellData<double> >(patch.getPatchData(d_rho, getDataContext())));

  tbox::for_all2<tbox::policy::parallel>(patch.getBox(), [=] __host__ __device__ (int k, int j) {
    rho(j,k) = 0.0;
  });
}

double
Stencil::computeStableDtOnPatch(
  hier::Patch& patch,
  const bool initial_time,
  const double dt_time)
{

  // calc dt

  return 0.5;
}

void
Stencil::computeFluxesOnPatch(
  hier::Patch& patch,
  const double time,
  const double dt)
{
  CellView<double, 2> rho(BOOST_CAST<pdat::CellData<double> >(patch.getPatchData(d_rho, getDataContext())));

  const int level = patch.getPatchLevelNumber();

  tbox::for_all2<tbox::policy::parallel>(patch.getBox(), [=] __host__ __device__ (int k, int j) {
    rho(j,k) = level;
  });
}

void
Stencil::conservativeDifferenceOnPatch(
  hier::Patch& patch,
  const double time,
  const double dt,
  bool at_syncronization)
{
  /* 
   * no-op
   */
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
}

void
Stencil::setPhysicalBoundaryConditions(
  hier::Patch& patch,
  const double fill_time,
  const hier::IntVector& ghost_width_to_fill)
{
}

void
Stencil::postprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{
}

 void
 Stencil::postprocessCoarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio)
{
}

 void
 Stencil::readDirichletBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
}

 void
 Stencil::readNeumannBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
}

void
Stencil::registerVisItDataWriter(boost::shared_ptr<appu::VisItDataWriter> viz_writer)
{
  d_visit_writer = viz_writer;
}
