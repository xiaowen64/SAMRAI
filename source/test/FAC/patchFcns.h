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
