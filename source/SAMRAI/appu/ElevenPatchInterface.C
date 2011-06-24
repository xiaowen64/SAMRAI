/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI interface to Eleven library 
 *
 ************************************************************************/
//
// THIS CLASS IS CURRENTLY EMPTY BECAUSE INTERFACES TO ELEVEN ARE
// STILL BEING DEVELOPED.  PLEASE CONTACT SAMRAI DEVELOPERS IF YOU
// ARE INTERESTED IN USING ELEVEN INTERFACES.
//
#ifndef included_ElevenPatchInterface_C
#define included_ElevenPatchInterface_C

#include "SAMRAI/appu/ElevenPatchInterface.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/appu/CutCell.h"
#include "SAMRAI/appu/EmbeddedBoundaryDefines.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/appu/ElevenPatchInterface.I"
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
ElevenPatchInterface::ElevenPatchInterface(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db)
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
ElevenPatchInterface::~ElevenPatchInterface()
{

}

/*
 *******************************************************************
 *
 *  Return eleven patch information
 *
 *******************************************************************
 */
void ElevenPatchInterface::calculateCutCellInfo(
   tbox::Pointer<hier::Patch>& patch,
   const int cell_flag_data_id,
   const int cell_vol_data_id,
   const int node_flag_data_id,
   const int cutcell_data_id)
{
   NULL_USE(patch);
   NULL_USE(cell_flag_data_id);
   NULL_USE(cell_vol_data_id);
   NULL_USE(node_flag_data_id);
   NULL_USE(cutcell_data_id);

   TBOX_ERROR(
      d_object_name << ":Unable to use ELEVEN."
      <<
      "\nPlease contact SAMRAI developers to find out more"
      << "\ninformation." << std::endl);
}

/*
 *******************************************************************
 *
 *  Check whether a particular node is inside the shape
 *
 *******************************************************************
 */
bool ElevenPatchInterface::isInside(
   const double* xyz) const
{
   NULL_USE(xyz);

   TBOX_ERROR(
      d_object_name << ":Unable to use ELEVEN."
      <<
      "\nPlease contact SAMRAI developers to find out more"
      << "\ninformation." << std::endl);

   return true;
}

/*
 *******************************************************************
 *
 *  Mark a SET of nodes on a patch as being inside or outside the
 *  shape.
 *
 *     nx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] is the dimensions of the patch
 *     dx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] is the delta x in each direction of the patch
 *     origin[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] is the location of the lower left node of the
 *                 patch
 *     inout[nx[0]*nx[1]*nx[2]] is the array holding whether the
 *                 nodes are inside or outside (1=inside, 0=outside)
 *
 *******************************************************************
 */
void ElevenPatchInterface::isInside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* inout) const
{
   NULL_USE(nx);
   NULL_USE(dx);
   NULL_USE(origin);
   NULL_USE(inout);

   TBOX_ERROR(
      d_object_name << ":Unable to use ELEVEN."
      <<
      "\nPlease contact SAMRAI developers to find out more"
      << "\ninformation." << std::endl);
}

/*
 *******************************************************************
 *
 *  Read info from input
 *
 *******************************************************************
 */
void ElevenPatchInterface::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
   NULL_USE(db);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   TBOX_ERROR(
      d_object_name << ":Unable to use ELEVEN."
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
void ElevenPatchInterface::printClassData(
   std::ostream& os) const
{
   os << "d_object_name = " << d_object_name << std::endl;
}

}
}
#endif
