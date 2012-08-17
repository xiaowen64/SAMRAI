/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating Cartesian 2d boundary data
 *
 ************************************************************************/

#ifndef included_appu_CartesianBoundaryUtilities2
#define included_appu_CartesianBoundaryUtilities2

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include "boost/shared_ptr.hpp"
#include <string>

namespace SAMRAI {
namespace appu {

/*!
 * @brief Class CartesianBoundaryUtilities2 is a utility class that
 * simplifies the implementation of simple physical boundary data in
 * 2 spatial dimensions.  It contains routines for reading boundary data
 * information from input files, applying those boundary conditions,
 * and error checking boundary data.  These routines apply to the
 * case of cell-centered double data only.  One may use all of these
 * capabilities, or use the input reading, boundary setting, and error
 * checking routines independently.
 *
 * <b> Input Parameters </b>
 *
 * To use the boundary condition input reading capabilities, the format
 * of the input file section containing the boundary information must
 * be as described next.  Boundary node and edge entries are only
 * required for those that are not filled automatically when periodic
 * conditions apply.
 *
 * The boundary condition for edge "*" is provided in a section as follows:
 *
 * @code
 *    boundary_edge_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 * @endcode
 *
 * Allowable edge identifiers (i.e., values for "*") are: <br>
 *       xlo, xhi, ylo, yhi <br>
 * Supported edge boundary_condition string values are: <br>
 *       "FLOW", "REFLECT", "DIRICHLET", "NEUMANN"
 *
 * The boundary condition for node "*" is provided in a section as follows:
 *
 * @code
 *    boundary_node_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 * @endcode
 *
 * Allowable node identifiers (i.e., values for "*") are: <br>
 *       xlo_ylo, xhi_ylo, xlo_yhi, xhi_yhi <br>
 * Supported node boundary_condition string values are: <br>
 *       "XFLOW", "YFLOW",
 *       "XREFLECT", "YREFLECT",
 *       "XDIRICHLET", "YDIRICHLET",
 *       "XNEUMANN", "YNEUMANN"
 *
 * Note that node conditions must be consistent with adjacent edge conditions.
 *
 * See the include file CartesianBoundaryDefines.h for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.  If you choose to use the input reading
 * capabilities only and write your own boundary condition routines in
 * FORTRAN, you should note that the integer constants for the various
 * boundary condition types and locations are automatically "stuffed" into
 * FORTRAN common blocks.  This avoids potential problems with
 * inconsistencies between C++ and FORTRAN usage.  Please see the
 * FORTRAN include file cartbdryparams2d.i for details.
 *
 * @see appu::BoundaryUtilityStrategy2
 */

struct CartesianBoundaryUtilities2 {
public:
   /*!
    * Function to read 2d boundary data from input database.
    * The integer boundary condition types are placed in the integer
    * arrays supplied by the caller (typically, the concrete
    * BoundaryUtilityStrategy object provided).  When DIRICHLET or
    * NEUMANN conditions are specified, control is passed to the
    * BoundaryUtilityStrategy to read the boundary state data specific to
    * the problem.
    *
    * Errors will be reported and the program will abort whenever necessary
    * boundary condition information is missing in the input database, or
    * when the data read in is either unknown or inconsistent.  The periodic
    * domain information is used to determine which boundary edge or
    * node entries are not required from input.  Error checking
    * requires that node boundary conditions are consistent
    * with those specified for the edges.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param bdry_strategy user-defined object that reads DIRICHLET or NEUMANN
    *                      conditions
    * @param input_db      input database containing all boundary data
    * @param edge_conds    array into which integer edge boundary condition types
    *                      are read
    * @param node_conds    array into which integer node boundary condition types
    *                      are read
    * @param periodic      integer vector specifying which coordinate
    *                      directions are periodic (e.g., value returned from
    *                      GridGeometry2::getPeriodicShift())
    */
   static void
   getFromInput(
      BoundaryUtilityStrategy* bdry_strategy,
      const boost::shared_ptr<tbox::Database>& input_db,
      tbox::Array<int>& edge_conds,
      tbox::Array<int>& node_conds,
      const hier::IntVector& periodic);

   /*!
    * Function to fill 2d edge boundary values for a patch.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param varname             String name of variable (for error reporting).
    * @param vardata             Cell-centered patch data object to fill.
    * @param patch               hier::Patch on which data object lives.
    * @param ghost_width_to_fill Width of ghost region to fill.
    * @param bdry_edge_conds     tbox::Array of boundary condition types for patch edges.
    * @param bdry_edge_values    tbox::Array of boundary values for patch edges.
    */
   static void
   fillEdgeBoundaryData(
      const std::string& varname,
      const boost::shared_ptr<pdat::CellData<double> >& vardata,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const tbox::Array<int>& bdry_edge_conds,
      const tbox::Array<double>& bdry_edge_values);

   /*!
    * Function to fill 2d node boundary values for a patch.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param varname             String name of variable (for error reporting).
    * @param vardata             Cell-centered patch data object to fill.
    * @param patch               hier::Patch on which data object lives.
    * @param ghost_width_to_fill Width of ghost region to fill.
    * @param bdry_node_conds     tbox::Array of boundary condition types for patch nodes.
    * @param bdry_edge_values    tbox::Array of boundary values for patch edges.
    */
   static void
   fillNodeBoundaryData(
      const std::string& varname,
      const boost::shared_ptr<pdat::CellData<double> >& vardata,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const tbox::Array<int>& bdry_node_conds,
      const tbox::Array<double>& bdry_edge_values);

   /*!
    * Function that returns the integer edge boundary location
    * corresponding to the given node location and node boundary
    * condition.
    *
    * If the node boundary condition type or node location are unknown,
    * or the boundary condition type is inconsistant with the node location
    * an error results.
    *
    * @return Integer edge location for node location and boundary condition type.
    *
    * @param node_loc   Integer location for node.
    * @param node_btype Integer boundary condition type for node.
    */
   static int
   getEdgeLocationForNodeBdry(
      int node_loc,
      int node_btype);

   /*!
    * Function to check 2d boundary data for a patch data quantity on
    * a patch after it is set.  A warning message will be sent to log
    * file for each bad boundary value that is found.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @return Integer number of bad boundary values found.
    *
    * @param varname       String name of variable (for error reporting).
    * @param patch         hier::Patch on which data object lives.
    * @param data_id       hier::Patch data index for data on patch.
    * @param depth         Depth index of patch data to check.
    * @param gcw_to_check  Width of ghost region to check.
    * @param bbox          Boundary box to check.
    * @param bcase         Boundary condition type for given edge or node.
    * @param bstate        Boundary value that applies in DIRICHLET or NEUMANN case.
    */
   static int
   checkBdryData(
      const std::string& varname,
      const hier::Patch& patch,
      int data_id,
      int depth,
      const hier::IntVector& gcw_to_check,
      const hier::BoundaryBox& bbox,
      int bcase,
      const double& bstate);

private:
   static bool s_fortran_constants_stuffed;

   static void
   read2dBdryEdges(
      BoundaryUtilityStrategy* bdry_strategy,
      const boost::shared_ptr<tbox::Database>& input_db,
      tbox::Array<int>& edge_conds,
      const hier::IntVector& periodic);

   static void
   read2dBdryNodes(
      const boost::shared_ptr<tbox::Database>& input_db,
      const tbox::Array<int>& edge_conds,
      tbox::Array<int>& node_conds,
      const hier::IntVector& periodic);

   static void
   get2dBdryDirectionCheckValues(
      int& idir,
      int& offsign,
      int btype,
      int bloc,
      int bcase);

   static void
   stuff2dBdryFortConst();
};

}
}
#endif
