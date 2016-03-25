/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	PoissonSineSolution class implementation
*/

#include "SAMRAI_config.h"

#include "MDA_Access.h"
#include "arrayConversion.h"
#include "patchFcns.h"
#include "PoissonSineSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "CartesianPatchGeometry.h"
#include "tbox/Array.h"

using namespace SAMRAI;

PoissonSineSolution::PoissonSineSolution()
  : d_linear_coef(0.0),
    d_exact()
{
  int i;
  for ( i=0; i<2*NDIM; ++i ) {
    d_neumann_location[i] = false;
  }
  return;
}

PoissonSineSolution::PoissonSineSolution(
  const string &object_name ,
  tbox::Database &database ,
  ostream *out_stream ,
  ostream *log_stream
)
  : d_linear_coef(0.0),
    d_exact()
{
  int i;
  for ( i=0; i<2*NDIM; ++i ) {
    d_neumann_location[i] = false;
  }
  setFromDatabase(database);
  return;
}

PoissonSineSolution::~PoissonSineSolution()
{
  return;
}

void PoissonSineSolution::setFromDatabase( tbox::Database &database )
{
  string istr = database.getStringWithDefault( "SinusoidFcnControl", "{}" );
  istringstream ist(istr);
  ist >> d_exact;
  if ( database.isBool("neumann_locations") ) {
    tbox::Array<bool> neumann_locations =
      database.getBoolArray("neumann_locations");
    if ( neumann_locations.getSize() > 2*NDIM ) {
      TBOX_ERROR("'neumann_locations' should have at most " << 2*NDIM
		 << " entries in " << NDIM << "D.\n");
    }
    int i;
    for ( i=0; i<neumann_locations.getSize(); ++i ) {
      d_neumann_location[i] = neumann_locations[i];
    }
  }
  d_linear_coef = database.getDoubleWithDefault( "linear_coef",
						 d_linear_coef );
  return;
}

void PoissonSineSolution::setNeumannLocation( int location_index, bool flag )
{
  assert( location_index < 2*NDIM );
  d_neumann_location[location_index] = flag;
  return;
}

void PoissonSineSolution::setPoissonSpecifications(
  solv::PoissonSpecifications &sps,
  int C_patch_data_id,
  int D_patch_data_id ) const
{
  sps.setDConstant(1.0);
  sps.setCConstant(d_linear_coef);
  return;
}

void PoissonSineSolution::setGridData (
  hier::Patch<NDIM> &patch ,
  pdat::SideData<NDIM,double> &diffcoef_data ,
  pdat::CellData<NDIM,double> &ccoef_data ,
  pdat::CellData<NDIM,double> &exact_data ,
  pdat::CellData<NDIM,double> &source_data )
{

  hier::Box<NDIM> pbox = patch.getBox();
  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
    = patch.getPatchGeometry();

  /* Set linear source coefficients */
  // ccoef_data.getArrayData().fill(d_linear_coef);

  /* Set diffusion coefficients. */
  // diffcoef_data.fillAll(1.0);

  /* Set source function and exact solution. */
  /*
    For the forcing function
    -( 1 + (nx^2 + ny^2) pi^2 ) * sin(pi*(nx*x+px))*sin(pi*(ny*y+py))
    the exact solution is sin(pi*(nx*x+px))*sin(pi*(ny*y+py))
  */
  setCellDataToSinusoid( exact_data,
			 patch,
			 d_exact );
  setCellDataToSinusoid( source_data,
			 patch,
			 d_exact );
  double npi[NDIM], ppi[NDIM];
  d_exact.getWaveNumbers(npi);
  d_exact.getPhaseAngles(ppi);
  const double source_scale = d_linear_coef -
#if NDIM == 2
		  ((npi[0]*npi[0]+npi[1]*npi[1])*M_PI*M_PI)
#endif
#if NDIM == 3
		  ((npi[0]*npi[0]+npi[1]*npi[1]+npi[2]*npi[2])*M_PI*M_PI)
#endif
    ;
  scaleArrayData( source_data.getArrayData() , source_scale );
  return;
}	// End patch loop.

ostream &operator<<( ostream &os, const PoissonSineSolution &r ) {
  os << r.d_exact << "\n";
  return os;
}

void PoissonSineSolution::setBcCoefs (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const tbox::Pointer< hier::Variable<NDIM> > &variable ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box ,
  const double fill_time ) const
{
  if ( bdry_box.getBoundaryType() != 1 ) {
    // Must be a face boundary.
    TBOX_ERROR("Bad boundary type in\n"
	       << "PoissonSineSolution::setBcCoefs \n");
  }

  const int location_index = bdry_box.getLocationIndex();

  // a is either 0 (Neumann) or 1 (Dirichlet).
  if ( !acoef_data.isNull() )
    acoef_data->fill( d_neumann_location[location_index] ? 0.0 : 1.0 , 0 );

  /*
    Get geometry information needed to compute coordinates
    of side centers.
  */
  hier::Box<NDIM> patch_box(patch.getBox());
  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
    = patch.getPatchGeometry();
  const double *xlo = patch_geom->getXLower();
  const double *xup = patch_geom->getXUpper();
  const double *dx = patch_geom->getDx();
  const hier::Box<NDIM> &box = bdry_box.getBox();
  hier::Index<NDIM> lower = box.lower();
  hier::Index<NDIM> upper = box.upper();

  if (gcoef_data) {
#if NDIM == 2
    hier::BoxIterator<NDIM> boxit( gcoef_data->getBox() );
    int i, j;
    double x, y;
    switch (location_index) {
    case 0:
      // min i edge
      x = xlo[0];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(1,0);
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  j = (*boxit)[1];
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  (*gcoef_data)(*boxit,0) = -slope(x,y);
	}
      }
      else {
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  j = (*boxit)[1];
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  (*gcoef_data)(*boxit,0) = d_exact(x,y);
	}
      }
      break;
    case 1:
      // max i edge
      x = xup[0];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(1,0);
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  j = (*boxit)[1];
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  (*gcoef_data)(*boxit,0) = slope(x,y);
	}
      }
      else {
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  j = (*boxit)[1];
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  (*gcoef_data)(*boxit,0) = d_exact(x,y);
	}
      }
      break;
    case 2:
      // min j edge
      y = xlo[1];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,1);
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  i = (*boxit)[0];
	  x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	  if(gcoef_data) (*gcoef_data)(*boxit,0) = -slope(x,y);
	}
      }
      else {
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  i = (*boxit)[0];
	  x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	  (*gcoef_data)(*boxit,0) = d_exact(x,y);
	}
      }
      break;
    case 3:
      // max j edge
      y = xup[1];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,1);
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  i = (*boxit)[0];
	  x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	  (*gcoef_data)(*boxit,0) = slope(x,y);
	}
      }
      else {
	if(gcoef_data) for ( ; boxit; boxit++ ) {
	  i = (*boxit)[0];
	  x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	  (*gcoef_data)(*boxit,0) = d_exact(x,y);
	}
      }
      break;
    default:
      TBOX_ERROR("Invalid location index in\n"
		 << "PoissonSineSolution<NDIM>::setBcCoefs");
    }
#endif

#if NDIM == 3
    MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > g_array;
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
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(1,0,0);
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( j=jbeg; j<=jend; ++j ) {
	    y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	    if(g_array) g_array(i,j,k) = -slope(x,y,z);
	  }
	}
      }
      else {
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( j=jbeg; j<=jend; ++j ) {
	    y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    case 1:
      // max i side
      jbeg = box.lower()[1]; jend = box.upper()[1];
      kbeg = box.lower()[2]; kend = box.upper()[2];
      i = box.upper()[0];
      x = xup[0];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(1,0,0);
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( j=jbeg; j<=jend; ++j ) {
	    y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	    if(g_array) g_array(i,j,k) = slope(x,y,z);
	  }
	}
      }
      else {
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( j=jbeg; j<=jend; ++j ) {
	    y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    case 2:
      // min j side
      ibeg = box.lower()[0]; iend = box.upper()[0];
      kbeg = box.lower()[2]; kend = box.upper()[2];
      j = box.lower()[1] + 1;
      y = xlo[1];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,1,0);
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = -slope(x,y,z);
	  }
	}
      }
      else {
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    case 3:
      // max j side
      ibeg = box.lower()[0]; iend = box.upper()[0];
      kbeg = box.lower()[2]; kend = box.upper()[2];
      j = box.upper()[1];
      y = xup[1];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,1,0);
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = slope(x,y,z);
	  }
	}
      }
      else {
	for ( k=kbeg; k<=kend; ++k ) {
	  z = xlo[2] + dx[2]*(k-patch_box.lower()[2]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    case 4:
      // min k side
      ibeg = box.lower()[0]; iend = box.upper()[0];
      jbeg = box.lower()[1]; jend = box.upper()[1];
      k = box.lower()[2] + 1;
      z = xlo[2];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,0,1);
	for ( j=jbeg; j<=jend; ++j ) {
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = -slope(x,y,z);
	  }
	}
      }
      else {
	for ( j=jbeg; j<=jend; ++j ) {
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    case 5:
      // max k side
      ibeg = box.lower()[0]; iend = box.upper()[0];
      jbeg = box.lower()[1]; jend = box.upper()[1];
      k = box.upper()[2];
      z = xup[2];
      if ( d_neumann_location[location_index] ) {
	SinusoidFcn slope = d_exact.differentiate(0,0,1);
	for ( j=jbeg; j<=jend; ++j ) {
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = slope(x,y,z);
	  }
	}
      }
      else {
	for ( j=jbeg; j<=jend; ++j ) {
	  y = xlo[1] + dx[1]*(j-patch_box.lower()[1]+0.5);
	  for ( i=ibeg; i<=iend; ++i ) {
	    x = xlo[0] + dx[0]*(i-patch_box.lower()[0]+0.5);
	    if(g_array) g_array(i,j,k) = d_exact(x,y,z);
	  }
	}
      }
      break;
    default:
      TBOX_ERROR("Invalid location index in\n"
		 << "PoissonSineSolution::setBcCoefs");
    }
#endif
  }

  return;
}

/*
***********************************************************************
* This class uses analytical boundary condition, so it can            *
* an unlimited number of extensions past the corner of a patch.       *
***********************************************************************
*/
hier::IntVector<NDIM> PoissonSineSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector<NDIM>(1000);
}
