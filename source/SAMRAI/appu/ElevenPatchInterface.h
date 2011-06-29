/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI interface to Eleven library 
 *
 ************************************************************************/

#ifndef included_ElevenPatchInterface
#define included_ElevenPatchInterface

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"

#ifdef HAVE_ELEVEN
#include "model.hh"
#endif

namespace SAMRAI {
namespace appu {

/*!
 * @brief This class provides an interface to the Eleven geometry library by
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
 * \endverbatim
 *
 * @see appu::EmbeddedBoundaryGeometry
 * @see appu::CutCell
 */

class ElevenPatchInterface
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
   ElevenPatchInterface(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db);

   /*!
    * The destructor does nothing.
    */
   virtual ~ElevenPatchInterface();

   /*!
    * This method computes information about the cut cells on a patch.
    */
   void
   calculateCutCellInfo(
      tbox::Pointer<hier::Patch>& patch,
      const int cell_flag_data_id,
      const int cell_vol_data_id,
      const int node_flag_data_id,
      const int cutcell_data_id);

   /*!
    * For a given point, specify whether that point is inside (true)
    * or outside (false) the set of shapes.
    *
    * @param xyz  double array[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying coordinates.
    */
   bool
   isInside(
      const double* xyz) const;

   /*!
    * For a SET of nodes on a patch, with patch dimensions (number of
    * nodes in each direction) nx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE], origin[tbox::Dimension::MAXIMUM_DIMENSION_VALUE], and grid spacing
    * dx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE], set the integer inout[nx[0]*nx[1]*nx[2]] array to specify
    * whether nodes are inside or outside.
    *
    * Node is:
    *     INSIDE  - set inout[ijk] = 1
    *     OUTSIDE - set inout[ijk] = 0
    *
    * @param nx integer array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying number of points in each dir
    * @param dx double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying spacing of points in each dir
    * @param origin double array [tbox::Dimension::MAXIMUM_DIMENSION_VALUE] specifying origin of lower corner
    * @param inout int array dimensioned the total number of points
    *        (i.e. nx[0]*nx[1]*nx[2]).  This is an OUTPUT quantity.
    */
   void
   isInside(
      const int* nx,
      const double* dx,
      const double* origin,
      int* inout) const;

   /*!
    * Dump data to supplied stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const;

private:
   /*
    * Read name, database_file, and geom tolerance information from input.
    * All inputs are optional.  The default behavior starts the interactive
    * interface to allow a user to build some geometry.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db);

   std::string d_object_name;

   /*
    * Eleven objects used to represent the geometry
    */
#ifdef HAVE_ELEVEN
//   mutable ELEVEN::PointList points;
//   mutable ELEVEN::CompoundCurveList curves;
   mutable ELEVEN::CompoundCurve::ProjectionInfo pi;
   mutable ELEVEN::ArrayVector3D x, xp;

   mutable ELEVEN::Model * d_eleven_model;
#endif
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/appu/ElevenPatchInterface.I"
#endif
#endif // included_ElevenPatchInterface
