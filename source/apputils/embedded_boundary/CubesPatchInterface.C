//
// File:        CubesPatchInterface.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 609 $
// Modified:    $Date: 2005-09-13 15:15:49 -0700 (Tue, 13 Sep 2005) $
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
#include "tbox/IEEE.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"


#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif


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
   const string& object_name,
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
              << "\ninformation." << endl);

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
   assert(!db.isNull());
#endif

  TBOX_ERROR(d_object_name << ":Unable to use CUBES due to license issues."
              << "\nPlease contact SAMRAI developers to find out more"
              << "\ninformation." << endl);


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
   ostream& os) const
{
   os << "d_object_name = " << d_object_name << endl;
}


}
}
#endif

