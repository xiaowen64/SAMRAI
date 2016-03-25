/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	PoissonGaussianDiffcoefSolution class declaration
*/

#ifndef included_PoissonGaussianDiffcoefSolution
#define included_PoissonGaussianDiffcoefSolution

#include <string>

#include "SinusoidFcn.h"
#include "GaussianFcn.h"

#include "tbox/Database.h"

/*
  SAMRAI classes
*/
#include "CellData.h"
#include "SideData.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"


using namespace SAMRAI;



#if NDIM == 2
#define _coords_ x,y
#define _coordsdef_ double x, double y
#elif NDIM == 3
#define _coords_ x,y,z
#define _coordsdef_ double x, double y, double z
#endif


/*!
  @brief Specialized class to provide Gaussian-diffcoef
  solution-specific stuff.

  The following variable coefficient problem is contrived
  to test the code on variable coefficients.

  The diffusion coefficients are
  @f[ D_x = D_y = e^{-\lambda |r-r_0|^2} @f]
  where @f$ r_0 @f$ is the center of the Gaussian.
  The exact solution
  @f[ u = sin(k_x x + \phi_x) sin(k_y y + \phi_y) sin(k_z z + \phi_z) @f]
  Source term (derived by substituting the exact solution into the PDE) is
  @f[
    e^{\lambda |r-r_0|^2}
    \left\{
      k^2 s_x s_y s_z
      - 2\lambda \left[ k_x(x-x_0) c_x s_y s_z
                      + k_y(y-y_0) s_x c_y s_z
                      + k_z(z-z_0) s_x s_y c_z
                 \right ]
      \right\}
  @f]
  where @f$ k^2 = k_x^2 + k_y^2 @f$,
  @f$ s_x = sin(k_x x + \phi_x) @f$,
  @f$ s_y = sin(k_y y + \phi_y) @f$,
  @f$ s_z = sin(k_z z + \phi_z) @f$,
  @f$ c_x = cos(k_x x + \phi_x) @f$,
  @f$ c_y = cos(k_y y + \phi_y) @f$ and
  @f$ c_z = cos(k_z z + \phi_z) @f$.
*/
class PoissonGaussianDiffcoefSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonGaussianDiffcoefSolution();

  PoissonGaussianDiffcoefSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonGaussianDiffcoefSolution();

  void setFromDatabase( tbox::Database &database );

  void setPoissonSpecifications(
    /*! Object to set */ solv::PoissonSpecifications &sps,
    /*! C id, if used */ int C_patch_data_id,
    /*! D id, if used */ int D_patch_data_id  ) const;

  /*!
    @brief Set parameters living on grid.

    Ignored data are: ccoef_data
    because it is constant.
  */
  void setGridData( hier::Patch<NDIM> &patch ,
		    pdat::SideData<NDIM,double> &diffcoef_data ,
		    pdat::CellData<NDIM,double> &ccoef_data ,
		    pdat::CellData<NDIM,double> &exact_data ,
		    pdat::CellData<NDIM,double> &source_data );

  void setBcCoefs (
    tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
    tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
    const tbox::Pointer< hier::Variable<NDIM> > &variable ,
    const hier::Patch<NDIM> &patch ,
    const hier::BoundaryBox<NDIM> &bdry_box ,
    double fill_time=0.0 ) const;

  hier::IntVector<NDIM> numberOfExtensionsFillable() const;

  //! Compute exact solution for a given coordinate.
  double exactFcn ( _coordsdef_ ) const ;
  //! Compute source for a given coordinate.
  double sourceFcn ( _coordsdef_ ) const ;
  //! Compute diffusion coefficient for a given coordinate.
  double diffcoefFcn ( _coordsdef_ ) const ;


  friend ostream &operator<<( ostream &os,
			      const PoissonGaussianDiffcoefSolution &r );

private:

  //! @brief Gaussian component of solution and source.
  GaussianFcn d_gcomp;
  //! @brief Sine-Sine component of solution and source.
  SinusoidFcn d_sscomp;
  //! @brief Cosine-Sine component of solution and source.
  SinusoidFcn d_cscomp;
  //! @brief Sine-Cosine component of solution and source.
  SinusoidFcn d_sccomp;
  //@{
  double d_lambda, d_k[NDIM], d_p[NDIM], d_k2;
  //@}

};


#endif	// included_PoissonGaussianDiffcoefSolution
