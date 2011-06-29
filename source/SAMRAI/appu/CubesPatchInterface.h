/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Cubes embedded boundary shape 
 *
 ************************************************************************/

#ifndef included_CubesPatchInterface
#define included_CubesPatchInterface

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Timer.h"

namespace SAMRAI {
namespace appu {

/**
 * This class provides an interface to the Cubes geometry library, the
 * embedded boundary grid generator in the Cart3D package from NASA Ames
 * (http://people.nas.nasa.gov/~aftosmis/cart3d/).
 *
 * Use of this class first requires that you link with the Cubes library.
 * The package may be obtained by filling out a license agreement at the
 * website listed above. Next, you must specify a 'CubesPatchInterface {...}"
 * input entry in the input for the EmbeddedBoundaryGeometry class.  In
 * this input you specify the name of the triangulated surface grid file
 * that defines the surface geometry.
 *
 * Once an instance of the class is created, SAMRAI will pass patch information
 * to the "calculateCutCellInfo()" class whenever an embedded boundary is
 * constructed on a new level.  Cubes computes the list of cut cells on the
 * patch and stores them as IndexData, templated on the CutCell class.
 *
 * Required input keys and data types:
 *
 *    - \b surface_tri_file
 *       string indicating the name of the triangulated surface file that
 *       defines the surface geometry.
 *
 *    - \b max_levels
 *       int value indicating the maximum number of levels used to construct
 *       the ADT tree in Cubes.  In general, this should be equal to the
 *       maximum number of levels in the hierarchy.  Because it is mainly used
 *       for storage, it may be larger, but should never be smaller than the
 *       number of hierarchy levels.
 *
 * Optional input keys:
 *
 *    - \b  verbose
 *       boolean turns on the 'verbose' flag in cubes.
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *   CubesPatchInterface {
 *      surface_tri_file = "slc-rev.tri"
 *      max_levels       = 5
 *      verbose          = TRUE
 *   }
 *
 * \endverbatim
 *
 * @see appu::EmbeddedBoundaryGeometry
 * @see appu::CutCell
 */
class CubesPatchInterface
{
public:
   /*!
    * The constructor initializes the Cubes library and optionally starts
    * an interactive interface
    *
    * @param object_name name of object of this class
    * @param input_db    the input database which contains the name of the
    *                    surface geometry.
    * @param grid_geom   grid geometry
    * @param nghosts     number of ghost cells
    *
    */
   CubesPatchInterface(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer<geom::CartesianGridGeometry> grid_geom,
      hier::IntVector nghosts);

   /*!
    * The destructor does nothing.
    */
   virtual ~CubesPatchInterface();

   /*!
    * Compute set of cut cells on a patch.  This method makes calls to
    * the cubes library, passing in the patch extents, and returning
    * the set of cut cells on the patch.
    *
    * @param patch pointer to a patch on a hierarchy level
    * @param cell_flag_data_id descriptor id for integer array holding cell
    *        flag. Stored on every cell of the patch, designates cell as flow,
    *        solid, or cut.
    * @param cell_vol_data_id descriptor id for double array holding volume
    *        fraction. Stored on every cell of the patch.
    * @param cutcell_data_id descriptor id for IndexData<CutCell> data.
    */
   void
   calculateCutCellInfo(
      tbox::Pointer<hier::Patch>& patch,
      const int cell_flag_data_id,
      const int cell_vol_data_id,
      const int cutcell_data_id);

   /*!
    * Set whether to record areas and normal.  Some calculations do not
    * require this information.  By default, it is set TRUE.
    */
   void
   setRecordAreasAndNormal(
      const bool record_an);

   /*!
    * Dump data to supplied stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*
    * Read name, command_file, and interactive start information from input.
    * All inputs are optional.  The default behavior starts the interactive
    * interface to allow a user to build some geometry.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db);

   std::string d_object_name;

   /*
    * Timers interspersed throughout the class.
    */
   static tbox::Pointer<tbox::Timer> t_cubes;
   static tbox::Pointer<tbox::Timer> t_set_cutcells;
   static tbox::Pointer<tbox::Timer> t_set_cutareas;

   /*
    * Whether or not to save the areas and normals.
    */
   bool d_record_areas_and_normal;

   /*
    * Cubes objects used to represent the geometry
    */

   /*
    * Triangulated surface file name.
    */
   std::string d_surf_tri_file;

   /*
    * Domain lower/upper bounds.
    */
   float d_domain_xlo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   float d_domain_xhi[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Number of vertices on coarsest level.
    */
   int d_domain_nx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Max number of levels for Cubes.  Eventually, we should couple this
    * with the gridding algorithm but we'll keep it as a separate entry
    * for now.
    */
   int d_max_levels;

   /*
    * Whether or not to use Cubes in "verbose" mode.
    */
   bool d_cubes_verbose;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/appu/CubesPatchInterface.I"
#endif
#endif // included_CubesPatchInterface
