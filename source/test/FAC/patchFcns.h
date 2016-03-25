/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/patchFcns.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Misc patch functions used in FAC solver tests.
 */

#include "SinusoidFcn.h"
#include "QuarticFcn.h"
#include "Patch.h"
#include "RefineSchedule.h"
#include "ArrayData.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "OutersideData.h"

using namespace SAMRAI;

void scaleArrayData(
  pdat::ArrayData<NDIM,double> &ad,
  double scale );

void setArrayDataToConstant(
  pdat::ArrayData<NDIM,double> &ad,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom,
  double value );

void setArrayDataTo(
  pdat::ArrayData<NDIM,double> &ad,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom );

void setCellDataToSinusoid(
  pdat::CellData<NDIM,double> &cd ,
  const hier::Patch<NDIM> &patch ,
  const SinusoidFcn &fcn );
/*!
  \brief Set pdat::ArrayData to Michael's exact solution.
*/
void setArrayDataToPerniceExact(
  pdat::ArrayData<NDIM,double> &ad,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom );
/*!
  \brief Set pdat::ArrayData to Michael's source function.
*/
void setArrayDataToPerniceSource(
  pdat::ArrayData<NDIM,double> &ad,
  const geom::CartesianPatchGeometry<NDIM> &patch_geom );

void setCellDataToQuartic(
  pdat::CellData<NDIM,double> &cd ,
  const hier::Patch<NDIM> &patch ,
  const QuarticFcn &fcn );
