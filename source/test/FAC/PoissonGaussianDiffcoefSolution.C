/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/PoissonGaussianDiffcoefSolution.C $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	PoissonGaussianDiffcoefSolution class implementation
*/

#include "SAMRAI_config.h"

#include "arrayConversion.h"
#include "patchFcns.h"
#include "PoissonGaussianDiffcoefSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "CartesianPatchGeometry.h"

using namespace SAMRAI;

PoissonGaussianDiffcoefSolution::PoissonGaussianDiffcoefSolution()
{
  return;
}

PoissonGaussianDiffcoefSolution::PoissonGaussianDiffcoefSolution(
  const string &object_name
, tbox::Database &database
, /*! Standard output stream */ ostream *out_stream
, /*! Log output stream */ ostream *log_stream
)
{
  setFromDatabase(database);
  return;
}

PoissonGaussianDiffcoefSolution::~PoissonGaussianDiffcoefSolution()
{
  return;
}

void PoissonGaussianDiffcoefSolution::setFromDatabase( tbox::Database &database )
{
  string istr;
  int i;
  {
    // Get the gaussian component.
    istr = database.getStringWithDefault( "GaussianFcnControl", "{}" );
    istringstream ist(istr);
    ist >> d_gcomp;
    d_lambda = d_gcomp.getLambda();
  }
  {
    // Get the sine-sine component.
    istr = database.getStringWithDefault( "SinusoidFcnControl", "{}" );
    istringstream ist(istr);
    ist >> d_sscomp;
    d_sscomp.getWaveNumbers(d_k);
    d_sscomp.getPhaseAngles(d_p);
    d_k2 = 0;
    for ( i=0; i<NDIM; ++i ) {
      d_k[i] *= M_PI;
      d_p[i] *= M_PI;
      d_k2 += d_k[i]*d_k[i];
    }
  }
  {
    // Compute the cosine-sine component.
    d_cscomp = d_sscomp;
    double new_phase_angles[NDIM];
    d_cscomp.getPhaseAngles(new_phase_angles);
    new_phase_angles[0] -= 0.5;
    d_cscomp.setPhaseAngles(new_phase_angles);
  }
  return;
}

double PoissonGaussianDiffcoefSolution::diffcoefFcn ( _coordsdef_
) const {
  return d_gcomp(_coords_);
}

double PoissonGaussianDiffcoefSolution::exactFcn ( _coordsdef_
) const {
  return d_sscomp(_coords_);
}

double PoissonGaussianDiffcoefSolution::sourceFcn ( _coordsdef_
) const {
  double rval;
  double trig_arg[NDIM];
  d_sscomp.getPhaseAngles(trig_arg);
  double gauss_ctr[NDIM];
  d_gcomp.getCenter(gauss_ctr);
#if NDIM == 2
  trig_arg[0] += d_k[0]*x;
  trig_arg[1] += d_k[1]*y;
  double sx=sin(trig_arg[0]), cx=cos(trig_arg[0]);
  double sy=sin(trig_arg[1]), cy=cos(trig_arg[1]);
  rval  = d_k[0] * (x-gauss_ctr[0]) * cx * sy;
  rval += d_k[1] * (y-gauss_ctr[1]) * sx * cy;
  rval *= 2*d_gcomp.getLambda();
  rval -= d_k2 * sx * sy;
  rval *= d_gcomp(_coords_);
#endif
#if NDIM == 3
  trig_arg[0] += d_k[0]*x;
  trig_arg[1] += d_k[1]*y;
  trig_arg[2] += d_k[2]*z;
  double sx=sin(trig_arg[0]), cx=cos(trig_arg[0]);
  double sy=sin(trig_arg[1]), cy=cos(trig_arg[1]);
  double sz=sin(trig_arg[2]), cz=cos(trig_arg[2]);
  rval  = d_k[0] * (x-gauss_ctr[0]) * cx * sy * sz;
  rval += d_k[1] * (y-gauss_ctr[1]) * sx * cy * sz;
  rval += d_k[2] * (z-gauss_ctr[2]) * sx * sy * cz;
  rval *= 2*d_gcomp.getLambda();
  rval -= d_k2 * sx * sy * sz;
  rval *= d_gcomp(_coords_);
#endif
  return rval;
}

void PoissonGaussianDiffcoefSolution::setPoissonSpecifications(
  solv::PoissonSpecifications &sps,
  int C_patch_data_id,
  int D_patch_data_id ) const
{
  sps.setDPatchDataId(D_patch_data_id);
  sps.setCZero();
  return;
}

void PoissonGaussianDiffcoefSolution::setGridData (
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
  int axis;
  {
    /* Set diffusion coefficients on each side at a time. */
    for ( axis=0; axis<NDIM; ++axis ) {
      double sl[NDIM]; // Like XLower, except for side.
      int j;
      for ( j=0; j<NDIM; ++j ) {
	sl[j] = j!=axis ? xl[j] + 0.5*h[j] : xl[j];
      }
      pdat::SideData<NDIM,double>::Iterator iter(patch.getBox(), axis);
      double _coords_;
      for ( ; iter; iter++ ) {
	const pdat::SideIndex<NDIM> &index = *iter;
	x = sl[0] + (index[0]-il[0])*h[0];
#if NDIM > 1
	y = sl[1] + (index[1]-il[1])*h[1];
#endif
#if NDIM > 2
	z = sl[2] + (index[2]-il[2])*h[2];
#endif
	diffcoef_data(index) = diffcoefFcn(_coords_);
      }
    }
  }
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
      y = sl[1] + (index[1]-il[1])*h[1];
#if NDIM > 2
      z = sl[2] + (index[2]-il[2])*h[2];
#endif
      exact_data(index) = exactFcn(_coords_);
      source_data(index) = sourceFcn(_coords_);
    }
  }

  return;
}	// End patch loop.

ostream &operator<<( ostream &os,
		     const PoissonGaussianDiffcoefSolution &r ) {
  os << r.d_gcomp << "\n";
  os << r.d_sscomp << "\n";
  os << r.d_cscomp << "\n";
  os << r.d_sccomp << "\n";
  return os;
}

void PoissonGaussianDiffcoefSolution::setBcCoefs (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &bcoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const tbox::Pointer< hier::Variable<NDIM> > &variable ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box ,
  const double fill_time ) const
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
	       << "PoissonGaussianDiffcoefSolution::setBcCoefs \n");
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
	       << "PoissonGaussianDiffcoefSolution::setBcCoefs");
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
	       << "PoissonGaussianDiffcoefSolution::setBcCoefs");
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
hier::IntVector<NDIM> PoissonGaussianDiffcoefSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector<NDIM>(1000);
}

#undef _coords_
