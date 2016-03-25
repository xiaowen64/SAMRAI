/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	PoissonPolynomialSolution class declaration
*/

#ifndef included_PoissonPolynomialSolution
#define included_PoissonPolynomialSolution


#include "QuarticFcn.h"


#include <string>
#include "tbox/Database.h"


/*
  SAMRAI classes
*/
#include "CellData.h"
#include "SideData.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"


using namespace SAMRAI;


/*!
  @brief Specialized class to provide polynomial-solution-specific stuff.
*/
class PoissonPolynomialSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonPolynomialSolution();

  PoissonPolynomialSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonPolynomialSolution();

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


  friend ostream &operator<<( ostream &os,
			      const PoissonPolynomialSolution &r );

private:

  QuarticFcn d_exact;
  QuarticFcn d_source;

};


#endif	// included_PoissonPolynomialSolution
