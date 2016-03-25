/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	PoissonMultigaussianSolution class implementation
*/

#include "SAMRAI_config.h"

#include "arrayConversion.h"
#include "patchFcns.h"
#include "PoissonMultigaussianSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "CartesianPatchGeometry.h"

#include "tbox/Utilities.h"

using namespace SAMRAI;

PoissonMultigaussianSolution::PoissonMultigaussianSolution()
#ifndef PACKAGE
  : d_ngauss(0)
#endif
{
  return;
}

PoissonMultigaussianSolution::PoissonMultigaussianSolution(
  const string &object_name
, tbox::Database &database
, /*! Standard output stream */ ostream *out_stream
, /*! Log output stream */ ostream *log_stream
)
#ifndef PACKAGE
  : d_ngauss(0)
#endif
{
  setFromDatabase(database);
  return;
}

PoissonMultigaussianSolution::~PoissonMultigaussianSolution()
{
  return;
}

void PoissonMultigaussianSolution::setFromDatabase( tbox::Database &database )
{
  string istr;
  {
    char singlegauss[100];
    sprintf( singlegauss, "GaussianFcnControl_%d", d_gauss_size );
    if ( ! database.isString(singlegauss) ) {
      TBOX_ERROR("You must have at least " << singlegauss << " defined in the\n"
		 "database for PoissonMultigaussianSolution.\n");
    }
    GaussianFcn gauss;
    do {
      istr = database.getString( singlegauss );
      istringstream ist(istr);
      ist >> gauss;
      d_gauss_append( gauss );
      sprintf( singlegauss, "GaussianFcnControl_%d", d_gauss_size );
    } while ( database.isString(singlegauss) );
  }
  return;
}

void PoissonMultigaussianSolution::setPoissonSpecifications(
  solv::PoissonSpecifications &sps,
  int C_patch_data_id,
  int D_patch_data_id ) const
{
  sps.setDConstant(1.0);
  sps.setCZero();
  return;
}

double PoissonMultigaussianSolution::exactFcn ( _coordsdef_
) const {
  double rval=0;
  d_gauss_const_iterator i;
  for ( i=d_gauss_begin; i!=d_gauss_end; ++i ) {
    rval += (*i)(_coords_);
  }
  return rval;
}

double PoissonMultigaussianSolution::sourceFcn ( _coordsdef_
) const {
  double rval=0;
  d_gauss_const_iterator i;
  for ( i=d_gauss_begin; i!=d_gauss_end; ++i ) {
    const GaussianFcn &gauss = *i;
    double gauss_ctr[NDIM];
    gauss.getCenter(gauss_ctr);
    double tval;
    tval = 4*gauss.getLambda()*( (x-gauss_ctr[0])*(x-gauss_ctr[0])
#if NDIM > 1
			       + (y-gauss_ctr[1])*(y-gauss_ctr[1])
#endif
#if NDIM > 2
			       + (z-gauss_ctr[2])*(z-gauss_ctr[2])
#endif
			       );
    tval += 2*NDIM;
    tval *= gauss(_coords_)*gauss.getLambda();
    rval += tval;
  }
  return rval;
}

void PoissonMultigaussianSolution::setGridData (
  hier::Patch<NDIM> &patch ,
  pdat::SideData<NDIM,double> &diffcoef_data ,
  pdat::CellData<NDIM,double> &ccoef_data ,
  pdat::CellData<NDIM,double> &exact_data ,
  pdat::CellData<NDIM,double> &source_data )
{

  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
    = patch.getPatchGeometry();

  const double *h = patch_geom->getDx();
  const double *xl = patch_geom->getXLower();
  const int *il = patch.getBox().lower();
  
  {
    /* Set cell-centered data. */
    double sl[NDIM]; // Like XLower, except for cell.
    int j;
    for ( j=0; j<NDIM; ++j ) {
      sl[j] = xl[j] + 0.5*h[j];
    }
    pdat::CellData<NDIM,double>::Iterator iter(patch.getBox());
    double _coords_;
    for ( ; iter; iter++ ) {
      const pdat::CellIndex<NDIM> &index = *iter;
      x = sl[0] + (index[0]-il[0])*h[0];
#if NDIM > 1
      y = sl[1] + (index[1]-il[1])*h[1];
#endif
#if NDIM > 2
      z = sl[2] + (index[2]-il[2])*h[2];
#endif
      exact_data(index) = exactFcn(_coords_);
      source_data(index) = sourceFcn(_coords_);
    }
  }

  return;
}	// End patch loop.

ostream &operator<<( ostream &os, const PoissonMultigaussianSolution &r ) {
  int i;
  os << "{\n";
  for ( i=0; i<r.d_gauss_size; ++i ) {
    os << "GaussianFcnControl_" << i << " " << r.d_gauss[i] << "\n";
  }
  os << "}\n";
  return os;
}

void PoissonMultigaussianSolution::setBcCoefs (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const tbox::Pointer< hier::Variable<NDIM> > &variable ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box ,
  const double fill_time ) const
{
  if ( acoef_data.isNull() && gcoef_data.isNull() ) {
    return;
  }

  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
    = patch.getPatchGeometry();
  /*
    Set to an inhomogeneous Dirichlet boundary condition.
  */
  hier::Box<NDIM> patch_box(patch.getBox());

  const double *xlo = patch_geom->getXLower();
  const double *xup = patch_geom->getXUpper();
  const double *dx = patch_geom->getDx();

  if ( bdry_box.getBoundaryType() != 1 ) {
    // Must be a face boundary.
    TBOX_ERROR("Bad boundary type in\n"
	       << "PoissonMultigaussianSolution::setBcCoefs \n");
  }
  const hier::Box<NDIM> &box = bdry_box.getBox();
  hier::Index<NDIM> lower = box.lower();
  hier::Index<NDIM> upper = box.upper();

#if NDIM == 2
  hier::BoxIterator<NDIM> boxit( acoef_data ?
			   acoef_data->getBox() : gcoef_data->getBox() );
  int i, j;
  double x, y;
  switch (bdry_box.getLocationIndex()) {
  case 0:
    // min i edge
    x = xlo[0];
    for ( ; boxit; boxit++ ) {
      j = (*boxit)[1];
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      if(acoef_data) (*acoef_data)(*boxit,0) = 1.0;
      if(gcoef_data) (*gcoef_data)(*boxit,0) = exactFcn(x,y);
    }
    break;
  case 1:
    // max i edge
    x = xup[0];
    for ( ; boxit; boxit++ ) {
      j = (*boxit)[1];
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      if(acoef_data) (*acoef_data)(*boxit,0) = 1.0;
      if(gcoef_data) (*gcoef_data)(*boxit,0) = exactFcn(x,y);
    }
    break;
  case 2:
    // min j edge
    y = xlo[1];
    for ( ; boxit; boxit++ ) {
      i = (*boxit)[0];
      x = xlo[1] + dx[1]*(i-patch_box.lower()[1]+0.5);
      if(acoef_data) (*acoef_data)(*boxit,0) = 1.0;
      if(gcoef_data) (*gcoef_data)(*boxit,0) = exactFcn(x,y);
    }
    break;
  case 3:
    // max j edge
    y = xup[1];
    for ( ; boxit; boxit++ ) {
      i = (*boxit)[0];
      x = xlo[1] + dx[1]*(i-patch_box.lower()[1]+0.5);
      if(acoef_data) (*acoef_data)(*boxit,0) = 1.0;
      if(gcoef_data) (*gcoef_data)(*boxit,0) = exactFcn(x,y);
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonMultigaussianSolution::setBcCoefs");
  }
#endif

#if NDIM == 3
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > a_array, g_array;
  if ( acoef_data ) a_array = arrayData2ArrayAccess( *acoef_data );
  if ( gcoef_data ) g_array = arrayData2ArrayAccess( *gcoef_data );
  int i, j, k, ibeg, iend, jbeg, jend, kbeg, kend;
  double x, y, z;
  switch (bdry_box.getLocationIndex()) {
  case 0:
    // min i side
    jbeg = box.lower()[1]; jend = box.upper()[1];
    kbeg = box.lower()[2]; kend = box.upper()[2];
    i = box.lower()[0] + 1;
    x = xlo[0];
    for ( k=kbeg; k<=kend; ++k ) {
      z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
      for ( j=jbeg; j<=jend; ++j ) {
	y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  case 1:
    // max i side
    jbeg = box.lower()[1]; jend = box.upper()[1];
    kbeg = box.lower()[2]; kend = box.upper()[2];
    i = box.upper()[0];
    x = xup[0];
    for ( k=kbeg; k<=kend; ++k ) {
      z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
      for ( j=jbeg; j<=jend; ++j ) {
	y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  case 2:
    // min j side
    ibeg = box.lower()[0]; iend = box.upper()[0];
    kbeg = box.lower()[2]; kend = box.upper()[2];
    j = box.lower()[1] + 1;
    y = xlo[1];
    for ( k=kbeg; k<=kend; ++k ) {
      z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
      for ( i=ibeg; i<=iend; ++i ) {
	x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  case 3:
    // max j side
    ibeg = box.lower()[0]; iend = box.upper()[0];
    kbeg = box.lower()[2]; kend = box.upper()[2];
    j = box.upper()[1];
    y = xup[1];
    for ( k=kbeg; k<=kend; ++k ) {
      z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
      for ( i=ibeg; i<=iend; ++i ) {
	x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  case 4:
    // min k side
    ibeg = box.lower()[0]; iend = box.upper()[0];
    jbeg = box.lower()[1]; jend = box.upper()[1];
    k = box.lower()[2] + 1;
    z = xlo[2];
    for ( j=jbeg; j<=jend; ++j ) {
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      for ( i=ibeg; i<=iend; ++i ) {
	x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  case 5:
    // max k side
    ibeg = box.lower()[0]; iend = box.upper()[0];
    jbeg = box.lower()[1]; jend = box.upper()[1];
    k = box.upper()[2];
    z = xup[2];
    for ( j=jbeg; j<=jend; ++j ) {
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      for ( i=ibeg; i<=iend; ++i ) {
	x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	if(a_array) a_array(i,j,k) = 1.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonMultigaussianSolution::setBcCoefs");
  }
#endif

  return;
}

/*
***********************************************************************
* This class uses analytical boundary condition, so it can            *
* an unlimited number of extensions past the corner of a patch.       *
***********************************************************************
*/
hier::IntVector<NDIM> PoissonMultigaussianSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector<NDIM>(1000);
}

#undef _coords_
