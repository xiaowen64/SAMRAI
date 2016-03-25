//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/embedded_boundary/CubesPatchInterface.C $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Cubes embedded boundary shape
//              
//
// THIS CLASS IS CURRENTLY EMPTY BECAUSE OF LICENSE ISSUES WITH CUBES.
// PLEASE CONTACT SAMRAI DEVELOPERS IF YOU ARE INTERESTED IN USING
// THE CUBES INTERFACES.
// 
#ifndef included_CubesPatchInterface_C
#define included_CubesPatchInterface_C

#include "CubesPatchInterface.h"

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CutCell.h"
#include "EmbeddedBoundaryDefines.h"
#include "IndexData.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"




namespace SAMRAI {
   namespace appu {

/*
*******************************************************************
*
*  Constructor
*
*******************************************************************
*/
template<int DIM> 
CubesPatchInterface<DIM>::CubesPatchInterface(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer<geom::CartesianGridGeometry<DIM> > grid_geom,
   hier::IntVector<DIM> nghosts)
{
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
template<int DIM> 
CubesPatchInterface<DIM>::~CubesPatchInterface<DIM>()
{

}


/*
*******************************************************************
*
*  Return cubes patch information
*
*******************************************************************
*/
template<int DIM> 
void CubesPatchInterface<DIM>::calculateCutCellInfo(
   tbox::Pointer<hier::Patch<DIM> >& patch,
   const int cell_flag_data_id,
   const int cell_vol_data_id,
   const int cutcell_data_id)
{
   TBOX_ERROR(d_object_name << ":Unable to use CUBES due to license issues."
              << "\nPlease contact SAMRAI developers to find out more"
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
template<int DIM>
void CubesPatchInterface<DIM>::setRecordAreasAndNormal(
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
template<int DIM> 
void CubesPatchInterface<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

  TBOX_ERROR(d_object_name << ":Unable to use CUBES due to license issues."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << std::endl);


}

/*
*******************************************************************
*
*  Dump (to supplied os) class information
*
*******************************************************************
*/
template<int DIM> 
void CubesPatchInterface<DIM>::printClassData(
   std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
}


}
}
#endif

