//
// File:        ElevenPatchInterface.h
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 605 $
// Modified:    $Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description: SAMRAI interface to Eleven library
//              
// 

#ifndef included_ElevenPatchInterface
#define included_ElevenPatchInterface

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif

#ifdef HAVE_ELEVEN
#include "model.hh"
using namespace ELEVEN;
#endif

namespace SAMRAI {
   namespace appu {


/**
 * This class provides an interface to the Eleven geometry library by
 * Kyle Chand in CASC.  The shapes over which the embedded boundary
 * is cut is defined through an XML database, the name of which is specified 
 * in the input file.  
 *
 * Use of this class first requires that you link with the Eleven library.
 * Next, you must specify a 'ElevenPatchInterface {...}" 
 * input entry in the input for the EmbeddedBoundaryGeometry class.  In
 * this input you specify the name of the XML database file that defines 
 * the surface geometry.  
 *
 * Once an instance of the class is created, SAMRAI will pass patch information
 * to the "calculateCutCellInfo()" class whenever an embedded boundary is
 * constructed on a new level.  The eleven interfaces may be used to 
 * compute the volume fractions through an inside/outside determination of
 * the nodes of the cell, and also may be used to compute the distance
 * of a particular node from the surface boundary. 
 *
 * Required input keys and data types:
 *
 *    - \b database_file       
 *       string indicating the name of the XML database file describing the
 *       geometry.
 * 
 *
 *
 * Optional input keys:
 * 
 *    - \b geom_tolerance
 *       double value indicating the error tolerance for determining
 *       geometric information relative to the cartesian grid.  For instance,
 *       a setting of 0.01 indicates the error tolerance will be 1% of the 
 *       supplied patch bounding box.  
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 * 
 *   ElevenPatchInterface {
 *      database_file = "circle_r.3.gzxml"
 *      geom_tolerance = 0.01 // 1% of bounding box
 *   }
 *
 * @see appu::EmbeddedBoundaryGeometry  
 * @see appu::CutCell  
 */
template<int DIM> class ElevenPatchInterface
{
public:
   /*!
    * The constructor initializes the eleven library and optionally starts 
    * an interactive interface
    *
    * @param object_name name of object of this class
    * @param input_db    the input database which contains radius and 
    *                    center specification.
    */
   ElevenPatchInterface(const string& object_name,
                        tbox::Pointer<tbox::Database> input_db);
   
   /*!
    * The destructor does nothing.
    */
   ~ElevenPatchInterface();

   /*!
    * This method computes information about the cut cells on a patch.
    */
   void calculateCutCellInfo(
      tbox::Pointer<hier::Patch<DIM> >& patch,
      const int cell_flag_data_id,
      const int cell_vol_data_id,
      const int node_flag_data_id,
      const int cutcell_data_id);
   
   /*!
    * For a given point, specify whether that point is inside (true)
    * or outside (false) the set of shapes.
    *
    * @param xyz  double array[DIM] specifying coordinates. 
    */
   bool isInside(const double* xyz) const;

   /*!
    * For a SET of nodes on a patch, with patch dimensions (number of
    * nodes in each direction) nx[DIM], origin[DIM], and grid spacing 
    * dx[DIM], set the integer inout[nx[0]*nx[1]*nx[2]] array to specify
    * whether nodes are inside or outside.  
    *
    * Node is:
    *     INSIDE  - set inout[ijk] = 1
    *     OUTSIDE - set inout[ijk] = 0
    *
    * @param nx integer array [DIM] specifying number of points in each dir
    * @param dx double array [DIM] specifying spacing of points in each dir
    * @param origin double array [DIM] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points 
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   void isInside(
      const int* nx,
      const double* dx,
      const double* origin,
      int* inout) const;
   
   
   /*!
    * Dump data to supplied stream.
    */
   void printClassData(ostream& os) const;

private:
   /*
    * Read name, database_file, and geom tolerance information from input.  
    * All inputs are optional.  The default behavior starts the interactive 
    * interface to allow a user to build some geometry.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   string d_object_name;

   /*
    * Eleven objects used to represent the geometry
    */
#ifdef HAVE_ELEVEN
//   mutable ELEVEN::PointList points;
//   mutable ELEVEN::CompoundCurveList curves;
   mutable ELEVEN::CompoundCurve::ProjectionInfo pi;
   mutable ELEVEN::ArrayVector3D x, xp;

   mutable ELEVEN::Model* d_eleven_model;
#endif

};   
 
}
}
#endif // included_ElevenPatchInterface


