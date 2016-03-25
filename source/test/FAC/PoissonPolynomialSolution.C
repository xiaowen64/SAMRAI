/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/test/FAC/PoissonPolynomialSolution.C $
  Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 2278 $
  Modified:	$LastChangedDate: 2008-07-07 14:51:53 -0700 (Mon, 07 Jul 2008) $
  Description:	PoissonPolynomialSolution class implementation
*/

#include "SAMRAI_config.h"

#include "MDA_Access.h"
#include "ArrayDataAccess.h"
#include "patchFcns.h"
#include "PoissonPolynomialSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "CartesianPatchGeometry.h"

using namespace SAMRAI;

PoissonPolynomialSolution::PoissonPolynomialSolution()
{
  return;
}

PoissonPolynomialSolution::PoissonPolynomialSolution(
  const string &object_name
, tbox::Database &database
, /*! Standard output stream */ ostream *out_stream
, /*! Log output stream */ ostream *log_stream
)
{
  setFromDatabase(database);
  return;
}

PoissonPolynomialSolution::~PoissonPolynomialSolution()
{
  return;
}

void PoissonPolynomialSolution::setFromDatabase( tbox::Database &database )
{
  string istr = database.getStringWithDefault( "QuarticFcnControl", "{}" );
  istringstream ist(istr);
  ist >> d_exact;
  /*
    Set the source to u_xx + u_yy + u_zz
  */
  d_source = d_exact.differentiate(2)
#if NDIM >= 2
           + d_exact.differentiate(0,2)
#endif
#if NDIM >= 3
           + d_exact.differentiate(0,0,2)
#endif
    ;

  return;
}

void PoissonPolynomialSolution::setPoissonSpecifications(
  solv::PoissonSpecifications &sps,
  int C_patch_data_id,
  int D_patch_data_id ) const
{
  sps.setDConstant(1.0);
  sps.setCZero();
  return;
}

void PoissonPolynomialSolution::setGridData (
  hier::Patch<NDIM> &patch ,
  pdat::SideData<NDIM,double> &diffcoef_data ,
  pdat::CellData<NDIM,double> &ccoef_data ,
  pdat::CellData<NDIM,double> &exact_data ,
  pdat::CellData<NDIM,double> &source_data )
{

  hier::Box<NDIM> pbox = patch.getBox();
  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
    = patch.getPatchGeometry();

  /* Set source function and exact solution. */
  setCellDataToQuartic( exact_data,
			patch,
			d_exact );
  setCellDataToQuartic( source_data,
			patch,
			d_source );
  return;
}	// End patch loop.

ostream &operator<<( ostream &os, const PoissonPolynomialSolution &r ) {
  os << r.d_exact << "\n";
  os << r.d_source << "\n";
  return os;
}

void PoissonPolynomialSolution::setBcCoefs (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &bcoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const tbox::Pointer< hier::Variable<NDIM> > &variable ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box,
  const double fill_time) const
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
	       << "PoissonPolynomialSolution::setBcCoefs \n");
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
      if(g_array) g_array[j-jbeg] = d_exact(x,y);
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
      if(g_array) g_array[j-jbeg] = d_exact(x,y);
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
      if(g_array) g_array[i-ibeg] = d_exact(x,y);
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
      if(g_array) g_array[i-ibeg] = d_exact(x,y);
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonPolynomialSolution::setBcCoefs");
  }
#endif

#if NDIM == 3
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > a_array, b_array, g_array;
  if ( acoef_data ) a_array = pdat::ArrayDataAccess::access( *acoef_data );
  if ( bcoef_data ) b_array = pdat::ArrayDataAccess::access( *bcoef_data );
  if ( gcoef_data ) g_array = pdat::ArrayDataAccess::access( *gcoef_data );
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
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
	if(g_array) g_array(i,j,k) = d_exact(x,y,z);
      }
    }
    break;
  default:
    TBOX_ERROR("Invalid location index in\n"
	       << "PoissonPolynomialSolution::setBcCoefs");
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
hier::IntVector<NDIM> PoissonPolynomialSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector<NDIM>(1000);
}
