//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/boundary/CartesianBoundaryUtilities2.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Utility routines for manipulating 2D Cartesian boundary data
//

#include "CartesianBoundaryUtilities2.h"

#include "CartesianBoundaryDefines.h"

#include "CartesianPatchGeometry.h"
#include "BoundaryBox.h"
#include "CellIndex.h"
#include "tbox/Utilities.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN 77 routines used in                 *
* boundary condition implementation.                                    *
*                                                                       *
*************************************************************************
*/

extern "C" {

void stufcartbdryloc2d_( const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const int& );

void stufcartbdrycond2d_( const int&, const int&, const int&,
                          const int&, const int&, const int&,
                          const int&, const int&, const int&,
                          const int&, const int&, const int& );

void getcartedgebdry2d_( const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&, 
                         const double*,
                         const int&,
                         const int&,
                         const double*,
                         double*,
                         const int& );

void getcartnodebdry2d_( const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const double*,
                         const int&,
                         const int&,
                         const double*,
                         double*,
                         const int& );

}

namespace SAMRAI {
    namespace appu {

bool CartesianBoundaryUtilities2::s_fortran_constants_stuffed = false;

/*
 * This function reads 2D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy 
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control 
 * is passed to the BoundaryStrategy to read the boundary state data 
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic 
 * domain information is used to determine which boundary edges or
 * node entries are not required from input.  Error checking requires 
 * that node boundary conditions are consistent with those 
 * specified for the edges.
 * 
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN data
 *    bdry_db .......... input database containing all boundary data
 *    edge_conds ....... array into which integer boundary conditions 
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions 
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate 
 *                       directions are periodic (value returned from
 *                       GridGeometry2::getPeriodicShift())
 */

void CartesianBoundaryUtilities2::readBoundaryInput(
   BoundaryUtilityStrategy* bdry_strategy,
   tbox::Pointer<tbox::Database> bdry_db, 
   tbox::Array<int>& edge_conds,
   tbox::Array<int>& node_conds,
   const hier::IntVector<2>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(bdry_strategy != (BoundaryUtilityStrategy*)NULL);
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(edge_conds.getSize() == NUM_2D_EDGES);
   TBOX_ASSERT(node_conds.getSize() == NUM_2D_NODES);
#endif

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   read2dBdryEdges(bdry_strategy,
                   bdry_db,
                   edge_conds,
                   periodic);

   read2dBdryNodes(bdry_db,
                   edge_conds,
                   node_conds,
                   periodic);

}

/*
 * Function to fill edge boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_edge_conds ...... array of boundary conditions for patch edges
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void CartesianBoundaryUtilities2::fillEdgeBoundaryData(
   const std::string& varname,
   tbox::Pointer< pdat::CellData<2,double> >& vardata,
   const hier::Patch<2>& patch,
   const hier::IntVector<2>& ghost_fill_width,
   const tbox::Array<int>& bdry_edge_conds,
   const tbox::Array<double>& bdry_edge_values)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(!vardata.isNull());
   TBOX_ASSERT(bdry_edge_conds.getSize() == NUM_2D_EDGES);
   TBOX_ASSERT(bdry_edge_values.getSize() == NUM_2D_EDGES*(vardata->getDepth()));
#endif

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom = 
      patch.getPatchGeometry();
   const double* dx = pgeom->getDx();

   const hier::Box<2>& interior = patch.getBox();
   const hier::Index<2>& ifirst = interior.lower();
   const hier::Index<2>& ilast  = interior.upper();

   const hier::IntVector<2>& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector<2> gcw_to_fill = hier::IntVector<2>::min(ghost_cells,
                                                      ghost_fill_width);

   const tbox::Array< hier::BoundaryBox<2> >& edge_bdry =
      pgeom->getCodimensionBoundary(EDGE2D_BDRY_TYPE);
   for (int i = 0; i < edge_bdry.getSize(); i++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(edge_bdry[i].getBoundaryType() == EDGE2D_BDRY_TYPE);
#endif

      int bedge_loc = edge_bdry[i].getLocationIndex();

      hier::Box<2> fill_box = pgeom->getBoundaryFillBox(edge_bdry[i],
                                                     interior,
                                                     gcw_to_fill);
      const hier::Index<2>& ibeg = fill_box.lower();
      const hier::Index<2>& iend = fill_box.upper();

      getcartedgebdry2d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         ghost_cells(0), ghost_cells(1), 
                         dx,
                         bedge_loc,
                         bdry_edge_conds[bedge_loc],
                         bdry_edge_values.getPointer(),
                         vardata->getPointer(),
                         vardata->getDepth());

   }

}

/*
 * Function to fill node boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void CartesianBoundaryUtilities2::fillNodeBoundaryData(
   const std::string& varname,
   tbox::Pointer< pdat::CellData<2,double> >& vardata,
   const hier::Patch<2>& patch,
   const hier::IntVector<2>& ghost_fill_width,
   const tbox::Array<int>& bdry_node_conds,
   const tbox::Array<double>& bdry_edge_values)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(!vardata.isNull());
   TBOX_ASSERT(bdry_node_conds.getSize() == NUM_2D_NODES);
   TBOX_ASSERT(bdry_edge_values.getSize() == NUM_2D_EDGES*(vardata->getDepth()));
#endif

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom = 
      patch.getPatchGeometry();
   const double* dx = pgeom->getDx();

   const hier::Box<2>& interior = patch.getBox();
   const hier::Index<2>& ifirst = interior.lower();
   const hier::Index<2>& ilast  = interior.upper();

   const hier::IntVector<2>& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector<2> gcw_to_fill = hier::IntVector<2>::min(ghost_cells,
                                                      ghost_fill_width);

   const tbox::Array< hier::BoundaryBox<2> >& node_bdry =
      pgeom->getCodimensionBoundary(NODE2D_BDRY_TYPE);
   for (int i = 0; i < node_bdry.getSize(); i++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(node_bdry[i].getBoundaryType() == NODE2D_BDRY_TYPE);
#endif

      int bnode_loc = node_bdry[i].getLocationIndex();

      hier::Box<2> fill_box = pgeom->getBoundaryFillBox(node_bdry[i],
                                                     interior,
                                                     gcw_to_fill);
      const hier::Index<2>& ibeg = fill_box.lower();
      const hier::Index<2>& iend = fill_box.upper();

      getcartnodebdry2d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         ghost_cells(0), ghost_cells(1),
                         dx,
                         bnode_loc,
                         bdry_node_conds[bnode_loc],
                         bdry_edge_values.getPointer(),
                         vardata->getPointer(),
                         vardata->getDepth());

   }

}

/*
 * Function that returns the integer edge boundary location 
 * corresponding to the given node location and node boundary
 * condition.  
 * 
 * If the node boundary condition type or node location are unknown, 
 * or the boundary condition type is inconsistant with the node location 
 * an error results.
 */

int CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
   int node_loc,
   int node_btype)
{

   int ret_edge = -1;

   switch (node_btype) {
      case XFLOW_BC:
      case XREFLECT_BC:
      case XDIRICHLET_BC:
      case XNEUMANN_BC:
      {
         if (node_loc == XLO_YLO_2D || node_loc == XLO_YHI_2D) { 
            ret_edge = XLO;
         } else {
            ret_edge = XHI;
         }
         break;
      }
      case YFLOW_BC:
      case YREFLECT_BC:
      case YDIRICHLET_BC:
      case YNEUMANN_BC:
      {
         if (node_loc == XLO_YLO_2D || node_loc == XHI_YLO_2D) {
            ret_edge = YLO;
         } else {
            ret_edge = YHI;
         }
         break;
      }
      default: {
         TBOX_ERROR("Unknown node boundary condition type = "
            << node_btype << " passed to \n"
            << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
            << std::endl); 
      }
   }

   if (ret_edge == -1) {
       TBOX_ERROR("Node boundary condition type = "
            << node_btype << " and node location = " << node_loc
            << "\n passed to "
            << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
            << " are inconsistant." << std::endl);
   }

   return(ret_edge);

}

/*
 * Function to check 2D boundary data filling.  Arguments are:
 *
 *    varname ..... name of variable (for error reporting)
 *    patch ....... patch on which boundary data to check lives
 *    data_id ..... patch data index on patch
 *    depth ....... depth index of data to check
 *    gcw_to_check. boundary ghost width to fill
 *    bbox ........ boundary box to check
 *    bcase ....... boundary condition case for edge or a node to check
 *    bstate ...... boundary state that applies when such a value is
 *                  required, such as when using Dirichlet conditions
 */

int CartesianBoundaryUtilities2::checkBdryData(
   const std::string& varname,
   const hier::Patch<2>& patch,
   int data_id,
   int depth,
   const hier::IntVector<2>& gcw_to_check,
   const hier::BoundaryBox<2>& bbox, 
   int bcase,
   double bstate)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(data_id >= 0);
   TBOX_ASSERT(depth >= 0);
#endif
  
   int num_bad_values = 0;

   int btype = bbox.getBoundaryType();
   int bloc  = bbox.getLocationIndex();

   tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom = 
      patch.getPatchGeometry();

   tbox::Pointer< pdat::CellData<2,double> > vardata = 
      patch.getPatchData(data_id);

   std::string bdry_type_str;
   if (btype == EDGE2D_BDRY_TYPE) {
      bdry_type_str = "EDGE";
   } else if (btype == NODE2D_BDRY_TYPE) {
      bdry_type_str = "NODE";
   } else {
      TBOX_ERROR("Unknown btype " << btype
         << " passed to CartesianBoundaryUtilities2::checkBdryData()! "
         << std::endl);
   }

   tbox::plog << "\n\nCHECKING 2D " << bdry_type_str << " BDRY DATA..." << std::endl;
   tbox::plog << "varname = " << varname << " : depth = " << depth << std::endl;
   tbox::plog << "bbox = " << bbox.getBox() << std::endl;
   tbox::plog << "btype, bloc, bcase = "
        << btype << ", = " << bloc << ", = " << bcase << std::endl;

   int idir;
   double valfact = 0.0, constval =0.0, dxfact=0.0;
   int offsign;

   get2dBdryDirectionCheckValues(idir, offsign,
                                 btype, bloc, bcase);

   if (btype == EDGE2D_BDRY_TYPE) {

      if (bcase == FLOW_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == REFLECT_BC) {
         valfact  = -1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == DIRICHLET_BC) {
         valfact  = 0.0;
         constval = bstate;
         dxfact   = 0.0;
      } else if (bcase == NEUMANN_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR("Unknown bcase " << bcase
            << " passed to CartesianBoundaryUtilities2::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   } else if (btype == NODE2D_BDRY_TYPE) {

      if (bcase == XFLOW_BC || bcase == YFLOW_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XREFLECT_BC || bcase == YREFLECT_BC) {
         valfact  = -1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XDIRICHLET_BC || bcase == YDIRICHLET_BC) {
         valfact  = 0.0;
         constval = bstate;
         dxfact   = 0.0;
      } else if (bcase == XNEUMANN_BC || bcase == YNEUMANN_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR("Unknown bcase " << bcase
            << " passed to CartesianBoundaryUtilities2::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   }

   hier::Box<2> gbox_to_check = 
      vardata->getGhostBox() * pgeom->getBoundaryFillBox(bbox,
                                                         patch.getBox(),
                                                         gcw_to_check);

   hier::Box<2> cbox = gbox_to_check;
   hier::Box<2> dbox = gbox_to_check;
   hier::Index<2> ifirst = vardata->getBox().lower();
   hier::Index<2> ilast  = vardata->getBox().upper();

   if (offsign == -1) {
      cbox.lower(idir) = ifirst(idir)-1;
      cbox.upper(idir) = ifirst(idir)-1;
      dbox.lower(idir) = ifirst(idir);
      dbox.upper(idir) = ifirst(idir);
   } else {
      cbox.lower(idir) = ilast(idir)+1;
      cbox.upper(idir) = ilast(idir)+1;
      dbox.lower(idir) = ilast(idir);
      dbox.upper(idir) = ilast(idir);
   }

   pdat::CellIterator<2> id(dbox);
   for (pdat::CellIterator<2> ic(cbox); ic; ic++) {
      double checkval = valfact * (*vardata)(id(), depth) + constval;
      pdat::CellIndex<2> check = ic();
      for (int p = 0; p < gbox_to_check.numberCells(idir); p++) {
         double offcheckval = checkval + dxfact * (p + 1);
         if ((*vardata)(check, depth) != offcheckval) {
            num_bad_values++;
            TBOX_WARNING("Bad " << bdry_type_str
                         << " boundary value for " << varname
                         << " found in cell " << check
                         << "\n   found = " << (*vardata)(check, depth)
                         << " : correct = " << offcheckval << std::endl);
         }
         check(idir) += offsign;
      }
      id++;
   }

   return(num_bad_values);

}

/*
 * Private function to read 2D edge boundary data from input database.
 */

void CartesianBoundaryUtilities2::read2dBdryEdges(
   BoundaryUtilityStrategy* bdry_strategy,
   tbox::Pointer<tbox::Database> bdry_db,
   tbox::Array<int>& edge_conds,
   const hier::IntVector<2>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(bdry_strategy != (BoundaryUtilityStrategy*)NULL);
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(edge_conds.getSize() == NUM_2D_EDGES);
#endif

   int num_per_dirs = 0;
   for (int id = 0; id < 2; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < 2) { // face boundary input required

      for (int s = 0; s < NUM_2D_EDGES; s++) {

         std::string bdry_loc_str;
         switch (s) {
            case XLO: { bdry_loc_str = "boundary_edge_xlo";
                        break; }
            case XHI: { bdry_loc_str = "boundary_edge_xhi";
                        break; }
            case YLO: { bdry_loc_str = "boundary_edge_ylo";
                        break; }
            case YHI: { bdry_loc_str = "boundary_edge_yhi";
                        break; }
            default: NULL_STATEMENT; 
         }

         bool need_data_read = true;
         if (num_per_dirs > 0) {
            if ( periodic(0) && (s == XLO || s == XHI) ) {
               need_data_read = false;
            } else if ( periodic(1) && (s == YLO || s == YHI) ) {
               need_data_read = false;
            }
         }

         if (need_data_read) {
            if (bdry_db->keyExists(bdry_loc_str)) {
               tbox::Pointer<tbox::Database> bdry_loc_db =
                  bdry_db->getDatabase(bdry_loc_str);
               if (!bdry_loc_db.isNull()) {
                  if (bdry_loc_db->keyExists("boundary_condition")) {
                     std::string bdry_cond_str =
                        bdry_loc_db->getString("boundary_condition");
                     if (bdry_cond_str == "FLOW") {
                        edge_conds[s] = FLOW_BC;
                     } else if (bdry_cond_str == "REFLECT") {
                        edge_conds[s] = REFLECT_BC;
                     } else if (bdry_cond_str == "DIRICHLET") {
                        edge_conds[s] = DIRICHLET_BC;
                        bdry_strategy->
                           readDirichletBoundaryDataEntry(bdry_loc_db,
                                                          bdry_loc_str,
                                                          s);
                     } else if (bdry_cond_str == "NEUMANN") {
                        edge_conds[s] = NEUMANN_BC;
                        bdry_strategy->
                           readNeumannBoundaryDataEntry(bdry_loc_db,
                                                        bdry_loc_str,
                                                        s);
                     } else {
                        TBOX_ERROR("Unknown edge boundary string = "
                           << bdry_cond_str << " found in input." << std::endl);
                     }
                  } else {
                     TBOX_ERROR("'boundary_condition' entry missing from "
                        << bdry_loc_str << " input database." << std::endl);
                  }
               }
            } else {
               TBOX_ERROR(bdry_loc_str
                  << " database entry not found in input." << std::endl);
            }
         } // if (need_data_read)

      } // for (int s = 0 ...

   } // if (num_per_dirs < 2)

}

/*
 * Private function to read 2D node boundary data from input database.
 */

void CartesianBoundaryUtilities2::read2dBdryNodes(
   tbox::Pointer<tbox::Database> bdry_db,
   const tbox::Array<int>& edge_conds,
   tbox::Array<int>& node_conds,
   const hier::IntVector<2>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(edge_conds.getSize() == NUM_2D_EDGES);
   TBOX_ASSERT(node_conds.getSize() == NUM_2D_NODES);
#endif

   int num_per_dirs = 0;
   for (int id = 0; id < 2; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < 1) { // node boundary data required 

      for (int s = 0; s < NUM_2D_NODES; s++) {

         std::string bdry_loc_str;
         switch (s) {
            case XLO_YLO_2D: { bdry_loc_str = "boundary_node_xlo_ylo";
                                break; }
            case XHI_YLO_2D: { bdry_loc_str = "boundary_node_xhi_ylo";
                                break; }
            case XLO_YHI_2D: { bdry_loc_str = "boundary_node_xlo_yhi";
                                break; }
            case XHI_YHI_2D: { bdry_loc_str = "boundary_node_xhi_yhi";
                                break; }
            default: NULL_STATEMENT;
         }

         if ( bdry_db->keyExists(bdry_loc_str) ) {
            tbox::Pointer<tbox::Database> bdry_loc_db = 
               bdry_db->getDatabase(bdry_loc_str);
            if (!bdry_loc_db.isNull()) {
               if (bdry_loc_db->keyExists("boundary_condition")) {
                  std::string bdry_cond_str =
                     bdry_loc_db->getString("boundary_condition");
                  if (bdry_cond_str == "XFLOW") {
                     node_conds[s] = XFLOW_BC;
                  } else if (bdry_cond_str == "YFLOW") {
                     node_conds[s] = YFLOW_BC;
                  } else if (bdry_cond_str == "XREFLECT") {
                     node_conds[s] = XREFLECT_BC;
                  } else if (bdry_cond_str == "YREFLECT") {
                     node_conds[s] = YREFLECT_BC;
                  } else if (bdry_cond_str == "XDIRICHLET") {
                     node_conds[s] = XDIRICHLET_BC;
                  } else if (bdry_cond_str == "YDIRICHLET") {
                     node_conds[s] = YDIRICHLET_BC;
                  } else if (bdry_cond_str == "XNEUMANN") {
                     node_conds[s] = XNEUMANN_BC;
                  } else if (bdry_cond_str == "YNEUMANN") {
                     node_conds[s] = YNEUMANN_BC;
                  } else {
                     TBOX_ERROR("Unknown node boundary string = "
                        << bdry_cond_str << " found in input." << std::endl);
                  }

                  std::string proper_edge;
                  std::string proper_edge_data;
                  bool no_edge_data_found = false;
                  if (bdry_cond_str == "XFLOW" ||
                      bdry_cond_str == "XDIRICHLET" ||
                      bdry_cond_str == "XNEUMANN" ||
                      bdry_cond_str == "XREFLECT") {
                     if (s == XLO_YLO_2D || s == XLO_YHI_2D) {
                        proper_edge = "XLO";
                        if (bdry_cond_str == "XFLOW" &&
                            edge_conds[XLO] != FLOW_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            edge_conds[XLO] != DIRICHLET_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            edge_conds[XLO] != NEUMANN_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            edge_conds[XLO] != REFLECT_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "REFLECT";
                        }
                     } else {
                        proper_edge = "XHI";
                        if (bdry_cond_str == "XFLOW" &&
                            edge_conds[XHI] != FLOW_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            edge_conds[XHI] != DIRICHLET_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            edge_conds[XHI] != NEUMANN_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            edge_conds[XHI] != REFLECT_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "REFLECT";
                        }
                     }
                  } else if (bdry_cond_str == "YFLOW" ||
                             bdry_cond_str == "YDIRICHLET" ||
                             bdry_cond_str == "YNEUMANN" ||
                             bdry_cond_str == "YREFLECT") {
                     if (s == XLO_YLO_2D || s == XHI_YLO_2D) {
                        proper_edge = "YLO";
                        if (bdry_cond_str == "YFLOW" &&
                            edge_conds[YLO] != FLOW_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            edge_conds[YLO] != DIRICHLET_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            edge_conds[YLO] != NEUMANN_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            edge_conds[YLO] != REFLECT_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "REFLECT";
                        }
                     } else {
                        proper_edge = "YHI";
                        if (bdry_cond_str == "YFLOW" &&
                            edge_conds[YHI] != FLOW_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            edge_conds[YHI] != DIRICHLET_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            edge_conds[YHI] != NEUMANN_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            edge_conds[YHI] != REFLECT_BC) {
                           no_edge_data_found = true;
                           proper_edge_data = "REFLECT";
                        }
                     }
                  } 
                  if (no_edge_data_found) {
                     TBOX_ERROR("Bdry condition " << bdry_cond_str
                        << " found for " << bdry_loc_str
                        << "\n but no " << proper_edge_data
                        << " data found for edge " << proper_edge << std::endl);
                  }

               } else {
                  TBOX_ERROR("'boundary_condition' entry missing from "
                     << bdry_loc_str << " input database." << std::endl);
               }
            }
         } else {
            TBOX_ERROR(bdry_loc_str
               << " database entry not found in input." << std::endl);
         }

      } // for (int s = 0 ...

   } // if (num_per_dirs < 1)

}

/*
 * Private function to get boundary orientation information for 
 * 2D boundary condition checking.  Called from checkBdryData().
 */

void CartesianBoundaryUtilities2::get2dBdryDirectionCheckValues(
   int& idir, 
   int& offsign, 
   int btype, 
   int bloc, 
   int bcase)
{

   std::string bdry_type_str;

   if (btype == EDGE2D_BDRY_TYPE) {

      bdry_type_str = "NODE";

      if (bloc == XLO || bloc == XHI) {
         idir = 0;
         if (bloc == XLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bloc == YLO || bloc == YHI) {
         idir = 1;
         if (bloc == YLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else {
         TBOX_ERROR("Unknown boundary location " << bloc
            << " passed to CartesianBoundaryUtilities2::checkBdryData()"
            << "\n for " << bdry_type_str << " boundary " << std::endl);
      }

   } else if (btype == NODE2D_BDRY_TYPE) {

      bdry_type_str = "NODE";

      if (bcase == XFLOW_BC || bcase == XREFLECT_BC ||
          bcase == XDIRICHLET_BC || bcase == XNEUMANN_BC) {
         idir = 0;
         if (bloc == XLO_YLO_2D || bloc == XLO_YHI_2D) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == YFLOW_BC || bcase == YREFLECT_BC ||
                 bcase == YDIRICHLET_BC || bcase == YNEUMANN_BC) {
         idir = 1;
         if (bloc == XLO_YLO_2D || bloc == XHI_YLO_2D) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } 

   } else {
      TBOX_ERROR("Unknown boundary type " << btype
         << " passed to CartesianBoundaryUtilities2::checkBdryData()"
         << "\n for " << bdry_type_str
         << " at location " << bloc << std::endl);
   }
  
}

/*
 * Private function to stuff 2D boundary contants into Fortran common blocks
 */

void CartesianBoundaryUtilities2::stuff2dBdryFortConst()
{
   stufcartbdryloc2d_(XLO, XHI, YLO, YHI,
                      XLO_YLO_2D, XHI_YLO_2D, XLO_YHI_2D, XHI_YHI_2D);
   stufcartbdrycond2d_(FLOW_BC,
                       XFLOW_BC, YFLOW_BC,
                       REFLECT_BC,
                       XREFLECT_BC, YREFLECT_BC,
                       DIRICHLET_BC,
                       XDIRICHLET_BC, YDIRICHLET_BC,
                       NEUMANN_BC,
                       XNEUMANN_BC, YNEUMANN_BC);
   s_fortran_constants_stuffed = true;
}

}
}
