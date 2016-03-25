/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/test/FAC/patchFcns.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2274 $
 * Modified:    $LastChangedDate: 2008-07-07 11:10:49 -0700 (Mon, 07 Jul 2008) $
 * Description: Misc patch functions used in FAC solver tests.
 */

#include "SAMRAI_config.h"

#include "MDA_Access.h"
#include "ArrayDataAccess.h"
#include "QuarticFcn.h"
#include "SinusoidFcn.h"

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
   MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = pdat::ArrayDataAccess::access(ad);
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
   MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = pdat::ArrayDataAccess::access(ad);
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
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = pdat::ArrayDataAccess::access(ad);
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
    t4 = pdat::ArrayDataAccess::access(cd.getArrayData());
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
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = pdat::ArrayDataAccess::access(ad);
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
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > t4 = pdat::ArrayDataAccess::access(ad);
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
    t4 = pdat::ArrayDataAccess::access(cd.getArrayData());
  setArrayDataToQuartic( t4 ,
			  cd.getGhostBox().lower() ,
			  cd.getGhostBox().upper() ,
			  cd.getBox().lower() ,
			  patch_geom->getXLower() ,
			  patch_geom->getDx() ,
			  fcn );
  return;
}
