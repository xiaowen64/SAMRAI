/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	AdaptivePoisson class declaration
*/

#ifndef included_AdaptivePoisson
#define included_AdaptivePoisson


#include "PoissonSpecifications.h"
#include "CellPoissonFACOps.h"
#include "PoissonSineSolution.h"
#include "PoissonPolynomialSolution.h"
#include "PoissonGaussianDiffcoefSolution.h"
#include "PoissonGaussianSolution.h"
#include "PoissonMultigaussianSolution.h"

#include <string>

#include "tbox/Pointer.h"
#include "tbox/Database.h"


/*
  SAMRAI classes
*/
#include "CartesianVizamraiDataWriter.h"
#include "VisItDataWriter.h"
#include "VisDerivedDataStrategy.h"
#include "CartesianCellDoubleConservativeLinearRefine.h"
#include "CartesianCellDoubleLinearRefine.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CartesianSideDoubleWeightedAverage.h"
#include "Box.h"
#include "CoarseFineBoundary.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "IntVector.h"
#include "StandardTagAndInitStrategy.h"
#include "CellDoubleConstantRefine.h"
#include "CellVariable.h"
#include "SideVariable.h"
#include "FaceVariable.h"
#include "NodeVariable.h"
#include "OutersideVariable.h"
#include "CartesianRobinBcHelper.h"
#include "FACPreconditioner.h"
#include "GhostCellRobinBcCoefs.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "RefinePatchStrategy.h"


using namespace SAMRAI;


/*!
  @brief Class to solve Poisson's equation on a SAMR grid.

  This class tests the FAC solver class solving
  Poisson's equation on a SAMR grid.

  This class inherits and implements virtual functions from
  mesh::StandardTagAndInitStrategy<NDIM> to initialize data
  on the SAMR grid.
*/
class AdaptivePoisson
  : public mesh::StandardTagAndInitStrategy<NDIM> ,
    public appu::VisDerivedDataStrategy<NDIM>
{

public:

  /*!
    @brief Constructor.

    Requirements:
    - the referenced objects

    Actions:
    - Set up private member data

    If you want standard output and logging,
    pass in valid pointers for those streams.
  */
  AdaptivePoisson(
    /*! Ojbect name */
    const string &object_name ,
    /*! Input database */
    tbox::Database &database ,
    /*! Standard output stream */ ostream *out_stream=NULL ,
    /*! Log output stream */ ostream *log_stream=NULL );



  //@{ @name mesh::StandardTagAndInitStrategy virtuals

public:

  /*!
    @brief Allocate and initialize data for a new level
    in the patch hierarchy.

    This is where you implement the code for initialize data on the
    grid.  Nevermind when it is called or where in the program that
    happens.  All the information you need to initialize the grid
    are in the arguments.

    @see mesh::StandardTagAndInitStrategy::initializeLevelData()
  */
  virtual void initializeLevelData (
    /*! Hierarchy to initialize */
    const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy ,
    /*! Level to initialize */
    const int level_number ,
    const double init_data_time ,
    const bool can_be_refined ,
    /*! Whether level is being introduced for the first time */
    const bool initial_time ,
    /*! Level to copy data from */
    const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level
      = tbox::Pointer< hier::BasePatchLevel<NDIM> >((0)) ,
    /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true );

  virtual void resetHierarchyConfiguration (
    /*! New hierarchy */
    tbox::Pointer<hier::BasePatchHierarchy<NDIM> > new_hierarchy ,
    /*! Coarsest level */ int coarsest_level ,
    /*! Finest level */ int finest_level );

  virtual void applyGradientDetector(
    const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation );

  //@}

#if 0
  //@{
  /*!
    @name Functions inherited from solv::CellPoissonFACOps<NDIM>
  */
  virtual void postprocessOneCycle(
    int iteration_num ,
    const solv::SAMRAIVectorReal<NDIM,double> &current_soln ,
    const solv::SAMRAIVectorReal<NDIM,double> &residual );
  //@}
#endif



  //@{ @name appu::VisDerivedDataStrategy<NDIM> virtuals

  virtual bool packDerivedDataIntoDoubleBuffer(
    double* buffer ,
    const hier::Patch<NDIM> &patch ,
    const hier::Box<NDIM> &region ,
    const string &variable_name ,
    int depth_id);

  //@}



public:

  /*!
    @brief Solve using FAC solver.

    Set up the linear algebra problem and use a
    solv::FACPreconditioner<NDIM> object to solve it.

    @param hierarchy the hierarchy to solve on
    @param max_cycle max number of FAC cycles to use
    @param pre_sweeps number of presmoothing sweeps to use
    @param pre_sweeps number of postsmoothing sweeps to use
    @param initial_u how to set the initial guess for u.
           A string is used so the option "random" can be
	   used.  If "random" is not used, set the string
	   to a floating point number.
  */
  int solvePoisson (
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy ,
    int max_cycles=10 ,
    double residual_tol=1e-6 ,
    int pre_sweeps=5 ,
    int post_sweeps=5 ,
    string initial_u=string("0.0") );

  /*!
    @brief Tell a plotter which data to write for this class.
  */
  int registerVariablesWithPlotter(
    appu::CartesianVizamraiDataWriter<NDIM> &viz_writer
    ) const;
  int registerVariablesWithPlotter(
    appu::VisItDataWriter<NDIM> &visit_writer
    );



  /*!
    @brief Compute the error of the current solution.

    Compute the @f$L_2@f$ and @f$L_\infty@f$ norms of the error,
    for each level and over all levels.
  */
  int computeError(
    /*! hierarchy */ const hier::PatchHierarchy<NDIM> &hierarchy ,
    /*! L2 norm */ double *l2norm ,
    /*! L-inf norm */ double *linorm ,
    /*! L2 norm on each level */ tbox::Array<double> &l2norms ,
    /*! L-inf norm on each level */ tbox::Array<double> &linorms ) const;


  /*!
    @brief Compute error estimator (for adaption or plotting).

    Computes in the box defined by @c estimate_data.
  */
  void computeAdaptionEstimate(
    pdat::CellData<NDIM,double> &estimate_data,
    const pdat::CellData<NDIM,double> &soln_cell_data) const;


private:
  string d_name;
  tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;

  //@{
  /*!
    @name Major algorithm objects.
  */
private:

  solv::CellPoissonFACOps<NDIM> d_fac_ops;

  solv::FACPreconditioner<NDIM> d_fac_preconditioner;

  //@}

  //@{
private:
  /*!
    @name Private state variables for solution.
  */

  /*!
    @brief Context for persistent data.
  */
  hier::VariableContext d_context_persistent;

  /*!
    @brief Context for scratch data.
  */
  hier::VariableContext d_context_scratch;

  /*!
    @brief Diffusion coefficient.
  */
  pdat::SideVariable<NDIM,double> d_diffcoef;

  /*!
    @brief Flux.
  */
  pdat::SideVariable<NDIM,double> d_flux;

  /*!
    @brief Scalar solution of Poisson's equation.
  */
  pdat::CellVariable<NDIM,double> d_scalar;

  /*!
    @brief Source for Poisson's equation.
  */
  pdat::CellVariable<NDIM,double> d_constant_source;

  /*!
    @brief Linear source operator for linear system.
  */
  pdat::CellVariable<NDIM,double> d_ccoef;

  /*!
    @brief Right hand side for linear system.
  */
  pdat::CellVariable<NDIM,double> d_rhs;

  /*!
    @brief Exact solution.
  */
  pdat::CellVariable<NDIM,double> d_exact;

  /*!
    @brief Residual.
  */
  tbox::Pointer<pdat::CellVariable<NDIM,double> > d_resid;

  /*!
    @brief Vector weights.

    For cells not covered by a finer cell, the weight
    is the volume.  For cells that are, the weight is zero.
    This is used in computing norms on the AMR grid.
  */
  pdat::CellVariable<NDIM,double> d_weight;


  /*!
    @brief Saved variable-context index.

    Because we refer to them often, variable-context indices are saved.
    They are initialized in the constructor and never change.
    Each index represents a variable-context pair in this class.
    Thus the indices are @em not independent state variables.
    If we had to set them, we need the associated variable and
    context (and the variable database which manages the mapping).
    See the hier::VariableDatabase<NDIM> class for more into.
  */
  int d_scalar_persistent, d_diffcoef_persistent ,
    d_constant_source_persistent, d_weight_persistent,
    d_exact_persistent , d_rhs_scratch, d_resid_scratch,
    d_flux_persistent, d_ccoef_persistent;

  //@}



  //@{
private:
  /*!
    @name Output streams.
  */
  /*!
    @brief Output stream pointer.

    If set to NULL, no output.
  */
  ostream *d_ostream;

  /*!
    @brief Log stream pointer.

    If set to NULL, no logging.
  */
  ostream *d_lstream;
  //@}

  //@{
  /*!
    @name Arbitrary function controls.
  */
  //! Wave number in half-cycles
  // double d_npi[NDIM];
  //! Phase shift in half-cycles
  // double d_ppi[NDIM];
  //@}


  //@{
  /*!
    @name Miscellaneous.
  */
  string d_problem_name;
  //! @brief Poisson equation specifications.
  solv::PoissonSpecifications d_sps;
  //! @brief Things specific to the sinusoid solution
  PoissonSineSolution d_sine_solution;
  //! @brief Things specific to the Gaussian solution
  PoissonGaussianSolution d_gaussian_solution;
  //! @brief Things specific to the multi-Gaussian solution
  PoissonMultigaussianSolution d_multigaussian_solution;
  //! @brief Things specific to the polynomial solution
  PoissonPolynomialSolution d_polynomial_solution;
  //! @brief Things specific to the Gaussian coefficient solution
  PoissonGaussianDiffcoefSolution d_gaussian_diffcoef_solution;
  /*!
    @brief Generic xfer::RefinePatchStrategy implementation for Robin bc.
  */
  solv::CartesianRobinBcHelper<NDIM> d_robin_refine_patch;
  /*!
    @brief Physical bc coefficient strategy selecting one of the solutions'.
  */
  solv::RobinBcCoefStrategy<NDIM> *d_physical_bc_coef;
  //@}

  double d_adaption_threshold;

  int d_finest_plot_level;

  //@{
private:
  /*!
    @name Objects to help debugging.
  */
  tbox::Pointer<appu::VisItDataWriter<NDIM> > d_visit_writer;
  int d_finest_dbg_plot_ln;
  //@}

};


#endif	// included_AdaptivePoisson
