/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	PoissonGaussianSolution class declaration
*/

#ifndef included_PoissonGaussianSolution
#define included_PoissonGaussianSolution

#include <string>

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
  @brief Specialized class to provide Gaussian
  solution-specific stuff.

  The exact solution is
  @f[ u = e^{-\lambda |r-r_0|^2} @f]
  where @f$ r_0 @f$ is the center of the Gaussian.

  The diffusion coefficients are 1.

  Plugging these into the Poisson equation, we get
  the following source function
  @f[ 2 \lambda e^{\lambda |r-r_0|^2} ( 3 + 2 \lambda |r-r0|^2 ) @f]
*/
class PoissonGaussianSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonGaussianSolution();

  PoissonGaussianSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonGaussianSolution();

  void setFromDatabase( tbox::Database &database );

  void setPoissonSpecifications(
    /*! Object to set */ solv::PoissonSpecifications &sps,
    /*! C id, if used */ int C_patch_data_id,
    /*! D id, if used */ int D_patch_data_id  ) const;

  /*!
    @brief Set parameters living on grid.

    Ignored data are: diffcoef_data and ccoef_data
    because they are constant.
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


  friend ostream &operator<<( ostream &os, const PoissonGaussianSolution &r );


private:

  //! @brief Gaussian component of solution and source.
  GaussianFcn d_gauss;

};


#endif	// included_PoissonGaussianSolution
