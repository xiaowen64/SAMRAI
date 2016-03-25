/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/PoissonMultigaussianSolution.h $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	PoissonMultigaussianSolution class declaration
*/

#ifndef included_PoissonMultigaussianSolution
#define included_PoissonMultigaussianSolution

#include <string>

#include "GaussianFcn.h"

/*
  Explanation of using PACKAGE macro to change code
  when compiling under SAMRAI:
  - SAMRAI does not define PACKAGE while this file's
    native package does.
  - We don't use the vector template under SAMRAI
    because implicit instantiation is disabled there.
*/

#ifdef PACKAGE
#include <vector>
#endif
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
  @brief Specialized class to provide solution-specific stuff
  for a solution with multiple Gaussians.

  A single Gaussian exact solution is
  @f[ u = e^{-\lambda |r-r_0|^2} @f]
  where @f$ r_0 @f$ is the center of the Gaussian.

  The diffusion coefficients are 1.

  Plugging these into the Poisson equation, we get
  the following source function
  @f[ 2 \lambda e^{\lambda |r-r_0|^2} ( 3 + 2 \lambda |r-r0|^2 ) @f]

  The multi-Gaussian is a sum of Gaussian solutions.
*/
class PoissonMultigaussianSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonMultigaussianSolution();

  PoissonMultigaussianSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonMultigaussianSolution();

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
      tbox::Pointer<pdat::ArrayData<NDIM,double> > &bcoef_data ,
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


  friend ostream &operator<<( ostream &os, const PoissonMultigaussianSolution &r );


private:

  //! @brief Gaussian component of solution and source.
#ifdef PACKAGE
  vector<GaussianFcn > d_gauss;
#define d_gauss_size d_gauss.size()
#define d_gauss_begin d_gauss.begin()
#define d_gauss_end d_gauss.end()
#define d_gauss_append(ITEM) d_gauss.insert( d_gauss.end(), ITEM )
#define d_gauss_const_iterator vector<GaussianFcn >::const_iterator
#else
  GaussianFcn d_gauss[20];
  int d_ngauss;
#define d_gauss_size d_ngauss
#define d_gauss_begin d_gauss
#define d_gauss_end (d_gauss+d_ngauss)
#define d_gauss_append(ITEM) { d_gauss[d_ngauss]=ITEM; d_ngauss++; }
#define d_gauss_const_iterator const GaussianFcn *
#endif

};


#endif	// included_PoissonMultigaussianSolution
