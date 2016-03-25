/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/test/FAC/PoissonGaussianSolution.C $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1846 $
  Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
  Description:	PoissonGaussianSolution class implementation
*/

#include "SAMRAI_config.h"

#include "arrayConversion.h"
#include "patchFcns.h"
#include "PoissonGaussianSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "CartesianPatchGeometry.h"

using namespace SAMRAI;

PoissonGaussianSolution::PoissonGaussianSolution()
{
  return;
}

PoissonGaussianSolution::PoissonGaussianSolution(
  const string &object_name
, tbox::Database &database
, /*! Standard output stream */ ostream *out_stream
, /*! Log output stream */ ostream *log_stream
)
{
  setFromDatabase(database);
  return;
}

PoissonGaussianSolution::~PoissonGaussianSolution()
{
  return;
}

void PoissonGaussianSolution::setFromDatabase( tbox::Database &database )
{
  string istr;
  {
    // Get the gaussian component.
    istr = database.getStringWithDefault( "GaussianFcnControl", "{}" );
    istringstream ist(istr);
    ist >> d_gauss;
  }
  return;
}

void PoissonGaussianSolution::setPoissonSpecifications(
  solv::PoissonSpecifications &sps,
  int C_patch_data_id,
  int D_patch_data_id ) const
{
  sps.setDConstant(1.0);
  sps.setCZero();
  return;
}

double PoissonGaussianSolution::exactFcn ( _coordsdef_
) const {
  return d_gauss(_coords_);
}

double PoissonGaussianSolution::sourceFcn ( _coordsdef_
) const {
  double gauss_ctr[NDIM];
  d_gauss.getCenter(gauss_ctr);
  double rval;
  rval = 4*d_gauss.getLambda()*( (x-gauss_ctr[0])*(x-gauss_ctr[0])
#if NDIM > 1
			       + (y-gauss_ctr[1])*(y-gauss_ctr[1])
#endif
#if NDIM > 2
			       + (z-gauss_ctr[2])*(z-gauss_ctr[2])
#endif
				);
  rval += 2*NDIM;
  rval *= d_gauss(_coords_)*d_gauss.getLambda();
  return rval;
}

void PoissonGaussianSolution::setGridData (
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

ostream &operator<<( ostream &os, const PoissonGaussianSolution &r ) {
  os << r.d_gauss << "\n";
  return os;
}

void PoissonGaussianSolution::setBcCoefs (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &bcoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const tbox::Pointer< hier::Variable<NDIM> > &variable ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box) const
{
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
	       << "PoissonGaussianSolution::setBcCoefs \n");
  }
  const hier::Box<NDIM> &box = bdry_box.getBox();
  hier::Index<NDIM> lower = box.lower();
  hier::Index<NDIM> upper = box.upper();

#if NDIM == 2
  double *a_array = acoef_data ? acoef_data->getPointer() : NULL;
  double *b_array = bcoef_data ? bcoef_data->getPointer() : NULL;
  double *g_array = gcoef_data ? gcoef_data->getPointer() : NULL;
  int i, j, ibeg, iend, jbeg, jend;
  double x, y;
  switch (bdry_box.getLocationIndex()) {
  case 0:
    // min i edge
    jbeg = box.lower()[1]; jend = box.upper()[1];
    x = xlo[0];
    for ( j=jbeg; j<=jend; ++j ) {
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      if(a_array) a_array[j-jbeg] = 1.0;
      if(b_array) b_array[j-jbeg] = 0.0;
      if(g_array) g_array[j-jbeg] = exactFcn(x,y);
    }
    break;
  case 1:
    // max i edge
    jbeg = box.lower()[1]; jend = box.upper()[1];
    x = xup[0];
    for ( j=jbeg; j<=jend; ++j ) {
      y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
      if(a_array) a_array[j-jbeg] = 1.0;
      if(b_array) b_array[j-jbeg] = 0.0;
      if(g_array) g_array[j-jbeg] = exactFcn(x,y);
    }
    break;
  case 2:
    // min j edge
    ibeg = box.lower()[0]; iend = box.upper()[0];
    y = xlo[1];
    for ( i=ibeg; i<=iend; ++i ) {
      x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
      if(a_array) a_array[i-ibeg] = 1.0;
      if(b_array) b_array[i-ibeg] = 0.0;
      if(g_array) g_array[i-ibeg] = exactFcn(x,y);
    }
    break;
  case 3:
    // max j edge
    ibeg = box.lower()[0]; iend = box.upper()[0];
    y = xup[1];
    for ( i=ibeg; i<=iend; ++i ) {
      x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
      if(a_array) a_array[i-ibeg] = 1.0;
      if(b_array) b_array[i-ibeg] = 0.0;
      if(g_array) g_array[i-ibeg] = exactFcn(x,y);
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonGaussianSolution::setBcCoefs");
  }
#endif

#if NDIM == 3
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > a_array, b_array, g_array;
  if ( acoef_data ) a_array = arrayData2ArrayAccess( *acoef_data );
  if ( bcoef_data ) b_array = arrayData2ArrayAccess( *bcoef_data );
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
	if(b_array) b_array(i,j,k) = 0.0;
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
	if(b_array) b_array(i,j,k) = 0.0;
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
	if(b_array) b_array(i,j,k) = 0.0;
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
	if(b_array) b_array(i,j,k) = 0.0;
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
	if(b_array) b_array(i,j,k) = 0.0;
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
	if(b_array) b_array(i,j,k) = 0.0;
	if(g_array) g_array(i,j,k) = exactFcn(x,y,z);
      }
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonGaussianSolution::setBcCoefs");
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
hier::IntVector<NDIM> PoissonGaussianSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector<NDIM>(1000);
}

#undef _coords_
