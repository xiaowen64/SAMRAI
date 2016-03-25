/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/PoissonSineSolution.h $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	PoissonSineSolution class declaration
*/

#ifndef included_PoissonSineSolution
#define included_PoissonSineSolution


#include "SinusoidFcn.h"


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
  @brief Specialized class to provide sine-solution-specific stuff.
*/
class PoissonSineSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonSineSolution();

  PoissonSineSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonSineSolution();

  void setFromDatabase( tbox::Database &database );

  void setNeumannLocation( int location_index, bool flag=true );

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


  friend ostream &operator<<( ostream &os, const PoissonSineSolution &r );

private:

  bool d_neumann_location[2*NDIM];
  double d_linear_coef;
  SinusoidFcn d_exact;

};


#endif	// included_PoissonSineSolution
