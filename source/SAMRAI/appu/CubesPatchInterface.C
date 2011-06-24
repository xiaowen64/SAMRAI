/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Cubes embedded boundary shape 
 *
 ************************************************************************/
//
// THIS CLASS IS CURRENTLY EMPTY BECAUSE OF LICENSE ISSUES WITH CUBES.
// PLEASE CONTACT SAMRAI DEVELOPERS IF YOU ARE INTERESTED IN USING
// THE CUBES INTERFACES.
//
#ifndef included_CubesPatchInterface_C
#define included_CubesPatchInterface_C

#include "SAMRAI/appu/CubesPatchInterface.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/appu/CutCell.h"
#include "SAMRAI/appu/EmbeddedBoundaryDefines.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/appu/CubesPatchInterface.I"
#endif

namespace SAMRAI {
namespace appu {

tbox::Pointer<tbox::Timer> CubesPatchInterface::t_cubes;
tbox::Pointer<tbox::Timer> CubesPatchInterface::t_set_cutcells;
tbox::Pointer<tbox::Timer> CubesPatchInterface::t_set_cutareas;

/*
 *******************************************************************
 *
 *  Constructor
 *
 *******************************************************************
 */
CubesPatchInterface::CubesPatchInterface(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer<geom::CartesianGridGeometry> grid_geom,
   hier::IntVector nghosts)
{
   NULL_USE(grid_geom);
   NULL_USE(nghosts);

   d_object_name = object_name;

   getFromInput(input_db);

}

/*
 *******************************************************************
 *
 *  Destructor
 *
 *******************************************************************
 */
CubesPatchInterface::~CubesPatchInterface()
{

}

/*
 *******************************************************************
 *
 *  Return cubes patch information
 *
 *******************************************************************
 */
void CubesPatchInterface::calculateCutCellInfo(
   tbox::Pointer<hier::Patch>& patch,
   const int cell_flag_data_id,
   const int cell_vol_data_id,
   const int cutcell_data_id)
{
   NULL_USE(patch);
   NULL_USE(cell_flag_data_id);
   NULL_USE(cell_vol_data_id);
   NULL_USE(cutcell_data_id);

   TBOX_ERROR(
      d_object_name << ":Unable to use CUBES due to license issues."
      <<
      "\nPlease contact SAMRAI developers to find out more"
      << "\ninformation." << std::endl);

}

/*
 *******************************************************************
 *
 * Set whether or not to record areas and normal.  By default, this
 * is turned on.
 *
 *******************************************************************
 */
void CubesPatchInterface::setRecordAreasAndNormal(
   const bool record_an)
{
   d_record_areas_and_normal = record_an;
}

/*
 *******************************************************************
 *
 *  Read info from input
 *
 *******************************************************************
 */
void CubesPatchInterface::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
   NULL_USE(db);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   TBOX_ERROR(
      d_object_name << ":Unable to use CUBES due to license issues."
      <<
      "\nPlease contact SAMRAI developers to find out more"
      << "\ninformation." << std::endl);

}

/*
 *******************************************************************
 *
 *  Dump (to supplied os) class information
 *
 *******************************************************************
 */
void CubesPatchInterface::printClassData(
   std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
}

}
}
#endif
