#include "SAMRAI_config.h"

#include "MDA_Access.h"
#include "QuarticFcn.h"
#include "SinusoidFcn.h"

#include "arrayConversion.h"
#include "setArrayData.h"

#include "Patch.h"
#include "tbox/Pointer.h"
#include "CartesianPatchGeometry.h"
#include "Box.h"
#include "CellData.h"
#include "SideData.h"
#include "OutersideData.h"

using namespace SAMRAI;

/*!
  \file
  \brief AMR-unaware functions to operate on a given single patch,
  to support FAC Poisson solve.
*/

/*!
  \brief Scale pdat::ArrayData.
*/
void scaleArrayData (
  pdat::ArrayData<NDIM,double> &ad ,
  double scale )
{
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = arrayData2ArrayAccess(ad);
  setArrayDataToScaled( t4 ,
			ad.getBox().lower() ,
			ad.getBox().upper() ,
			scale );
  return;
}

/*!
  \brief Set pdat::ArrayData to a constant.
*/
void setArrayDataToConstant(
  pdat::ArrayData<NDIM,double> &ad ,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom ,
  double value )
{
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = arrayData2ArrayAccess(ad);
  setArrayDataToConstant( t4 ,
			  ad.getBox().lower() ,
			  ad.getBox().upper() ,
			  patch_geom.getXLower() ,
			  patch_geom.getXUpper() ,
			  patch_geom.getDx() ,
			  value );
  return;
}

/*!
  \brief Set pdat::ArrayData to the x coordinate.
*/
void setArrayDataTo(
  pdat::ArrayData<NDIM,double> &ad ,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom )
{
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = arrayData2ArrayAccess(ad);
  setArrayDataTo( t4 ,
		   ad.getBox().lower() ,
		   ad.getBox().upper() ,
		   patch_geom.getXLower() ,
		   patch_geom.getXUpper() ,
		   patch_geom.getDx() );
  return;
}

/*!
  \brief Set pdat::CellData<NDIM> to a sinusoid function.
*/
void setCellDataToSinusoid(
  pdat::CellData<NDIM,double> &cd ,
  const hier::Patch<NDIM> &patch ,
  const SinusoidFcn &fcn )
{
  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >
    patch_geom = patch.getPatchGeometry();
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >
    t4 = arrayData2ArrayAccess(cd.getArrayData());
  setArrayDataToSinusoid( t4 ,
			  cd.getGhostBox().lower() ,
			  cd.getGhostBox().upper() ,
			  cd.getBox().lower() ,
			  patch_geom->getXLower() ,
			  patch_geom->getDx() ,
			  fcn );
  return;
}

/*!
  \brief Set pdat::ArrayData to Michael's exact solution.
*/
void setArrayDataToPerniceExact(
  pdat::ArrayData<NDIM,double> &ad ,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom )
{
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = arrayData2ArrayAccess(ad);
  setArrayDataToPerniceExact( t4 ,
			      ad.getBox().lower() ,
			      ad.getBox().upper() ,
			      patch_geom.getXLower() ,
			      patch_geom.getXUpper() ,
			      patch_geom.getDx() );
  return;
}

/*!
  \brief Set pdat::ArrayData to Michael's source function.
*/
void setArrayDataToPerniceSource(
  pdat::ArrayData<NDIM,double> &ad ,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom )
{
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = arrayData2ArrayAccess(ad);
  setArrayDataToPerniceSource( t4 ,
			       ad.getBox().lower() ,
			       ad.getBox().upper() ,
			       patch_geom.getXLower() ,
			       patch_geom.getXUpper() ,
			       patch_geom.getDx() );
  return;
}

/*!
  \brief Set pdat::ArrayData to a quartic function.
*/
void setCellDataToQuartic(
  pdat::CellData<NDIM,double> &cd ,
  const hier::Patch<NDIM> &patch ,
  const QuarticFcn &fcn )
{
  tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >
    patch_geom = patch.getPatchGeometry();
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >
    t4 = arrayData2ArrayAccess(cd.getArrayData());
  setArrayDataToQuartic( t4 ,
			  cd.getGhostBox().lower() ,
			  cd.getGhostBox().upper() ,
			  cd.getBox().lower() ,
			  patch_geom->getXLower() ,
			  patch_geom->getDx() ,
			  fcn );
  return;
}
