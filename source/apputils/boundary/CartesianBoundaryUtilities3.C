//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/boundary/CartesianBoundaryUtilities3.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Utility routines for manipulating 3D Cartesian boundary data
//

#include "CartesianBoundaryUtilities3.h"

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

void stufcartbdryloc3d_( const int&, const int&, 
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const int& );

void stufcartbdrycond3d_( const int&, 
                          const int&, const int&, const int&,
                          const int&, 
                          const int&, const int&, const int&,
                          const int&,
                          const int&, const int&, const int&,
                          const int&,
                          const int&, const int&, const int& );

void getcartfacebdry3d_( const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&, const int&,
                         const double*,
                         const int&,
                         const int&,
                         const double*,
                         double*, 
                         const int& );

void getcartedgebdry3d_( const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&, const int&,
                         const double*,
                         const int&,
                         const int&,
                         const double*,
                         double*,
                         const int& );

void getcartnodebdry3d_( const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&, const int&,
                         const double*,
                         const int&,
                         const int&,
                         const double*,
                         double*,
                         const int& );

}

namespace SAMRAI {
    namespace appu {

bool CartesianBoundaryUtilities3::s_fortran_constants_stuffed = false;

/*
 * This function reads 3D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy 
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control 
 * is passed to the BoundaryStrategy to read the boundary state data 
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic 
 * domain information is used to determine which boundary face, edge, or
 * node entries are not required from input.  Error checking requires 
 * that node and edge boundary conditions are consistent with those 
 * specified for the faces.
 * 
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN conditions
 *    bdry_db .......... input database containing all boundary data
 *    face_conds ....... array into which integer boundary conditions 
 *                       for faces are read
 *    edge_conds ....... array into which integer boundary conditions 
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions 
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate 
 *                       directions are periodic (value returned from
 *                       GridGeometry3::getPeriodicShift())
 */

void CartesianBoundaryUtilities3::readBoundaryInput(
   BoundaryUtilityStrategy* bdry_strategy,
   tbox::Pointer<tbox::Database> bdry_db, 
   tbox::Array<int>& face_conds,
   tbox::Array<int>& edge_conds,
   tbox::Array<int>& node_conds,
   const hier::IntVector<3>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(bdry_strategy != (BoundaryUtilityStrategy*)NULL);
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(face_conds.getSize() == NUM_3D_FACES);
   TBOX_ASSERT(edge_conds.getSize() == NUM_3D_EDGES);
   TBOX_ASSERT(node_conds.getSize() == NUM_3D_NODES);
#endif

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   read3dBdryFaces(bdry_strategy,
                   bdry_db,
                   face_conds,
                   periodic);

   read3dBdryEdges(bdry_db,
                   face_conds,
                   edge_conds,
                   periodic); 

   read3dBdryNodes(bdry_db,
                   face_conds,
                   node_conds,
                   periodic);

}

/*
 * Function to fill face boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_face_conds ...... array of boundary conditions for patch faces
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void CartesianBoundaryUtilities3::fillFaceBoundaryData(
   const std::string& varname,
   tbox::Pointer< pdat::CellData<3,double> >& vardata,
   const hier::Patch<3>& patch,
   const hier::IntVector<3>& ghost_fill_width,
   const tbox::Array<int>& bdry_face_conds,
   const tbox::Array<double>& bdry_face_values)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(!vardata.isNull());
   TBOX_ASSERT(bdry_face_conds.getSize() == NUM_3D_FACES);
   TBOX_ASSERT(bdry_face_values.getSize() == NUM_3D_FACES*(vardata->getDepth()));
#endif

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<3> > pgeom = 
      patch.getPatchGeometry();
   const double* dx = pgeom->getDx();

   const hier::Box<3>& interior = patch.getBox();
   const hier::Index<3>& ifirst = interior.lower();
   const hier::Index<3>& ilast  = interior.upper();

   const hier::IntVector<3>& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector<3> gcw_to_fill = hier::IntVector<3>::min(ghost_cells,
                                                      ghost_fill_width);

   const tbox::Array<hier::BoundaryBox<3> >& face_bdry =
      pgeom->getCodimensionBoundary(FACE3D_BDRY_TYPE);
   for (int i = 0; i < face_bdry.getSize(); i++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(face_bdry[i].getBoundaryType() == FACE3D_BDRY_TYPE);
#endif

      int bface_loc = face_bdry[i].getLocationIndex();

      hier::Box<3> fill_box = pgeom->getBoundaryFillBox(face_bdry[i],
                                                     interior,
                                                     gcw_to_fill);
      const hier::Index<3>& ibeg = fill_box.lower();
      const hier::Index<3>& iend = fill_box.upper();

      getcartfacebdry3d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ifirst(2), ilast(2),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         ibeg(2), iend(2),
                         ghost_cells(0), ghost_cells(1), ghost_cells(2),
                         dx,
                         bface_loc,
                         bdry_face_conds[bface_loc],
                         bdry_face_values.getPointer(),
                         vardata->getPointer(),
                         vardata->getDepth());

   }

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
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void CartesianBoundaryUtilities3::fillEdgeBoundaryData(
   const std::string& varname,
   tbox::Pointer< pdat::CellData<3,double> >& vardata,
   const hier::Patch<3>& patch,
   const hier::IntVector<3>& ghost_fill_width,
   const tbox::Array<int>& bdry_edge_conds,
   const tbox::Array<double>& bdry_face_values)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(!vardata.isNull());
   TBOX_ASSERT(bdry_edge_conds.getSize() == NUM_3D_EDGES);
   TBOX_ASSERT(bdry_face_values.getSize() == NUM_3D_FACES*(vardata->getDepth()));
#endif

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<3> > pgeom = 
      patch.getPatchGeometry();
   const double* dx = pgeom->getDx();

   const hier::Box<3>& interior = patch.getBox();
   const hier::Index<3>& ifirst = interior.lower();
   const hier::Index<3>& ilast  = interior.upper();

   const hier::IntVector<3>& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector<3> gcw_to_fill = hier::IntVector<3>::min(ghost_cells,
                                                      ghost_fill_width);

   const tbox::Array<hier::BoundaryBox<3> >& edge_bdry =
      pgeom->getCodimensionBoundary(EDGE3D_BDRY_TYPE);
   for (int i = 0; i < edge_bdry.getSize(); i++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(edge_bdry[i].getBoundaryType() == EDGE3D_BDRY_TYPE);
#endif

      int bedge_loc = edge_bdry[i].getLocationIndex();

      hier::Box<3> fill_box = pgeom->getBoundaryFillBox(edge_bdry[i],
                                                     interior,
                                                     gcw_to_fill);
      const hier::Index<3>& ibeg = fill_box.lower();
      const hier::Index<3>& iend = fill_box.upper();

      getcartedgebdry3d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ifirst(2), ilast(2),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         ibeg(2), iend(2),
                         ghost_cells(0), ghost_cells(1), ghost_cells(2),
                         dx,
                         bedge_loc,
                         bdry_edge_conds[bedge_loc],
                         bdry_face_values.getPointer(),
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
 *    bdry_face_values ..... array of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void CartesianBoundaryUtilities3::fillNodeBoundaryData(
   const std::string& varname,
   tbox::Pointer< pdat::CellData<3,double> >& vardata,
   const hier::Patch<3>& patch,
   const hier::IntVector<3>& ghost_fill_width,
   const tbox::Array<int>& bdry_node_conds,
   const tbox::Array<double>& bdry_face_values)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(!vardata.isNull());
   TBOX_ASSERT(bdry_node_conds.getSize() == NUM_3D_NODES);
   TBOX_ASSERT(bdry_face_values.getSize() == NUM_3D_FACES*(vardata->getDepth()));
#endif

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<3> > pgeom = 
      patch.getPatchGeometry();
   const double* dx = pgeom->getDx();

   const hier::Box<3>& interior = patch.getBox();
   const hier::Index<3>& ifirst = interior.lower();
   const hier::Index<3>& ilast  = interior.upper();

   const hier::IntVector<3>& ghost_cells = vardata->getGhostCellWidth();
 
   hier::IntVector<3> gcw_to_fill = hier::IntVector<3>::min(ghost_cells,
                                                      ghost_fill_width);

   const tbox::Array<hier::BoundaryBox<3> >& node_bdry =
      pgeom->getCodimensionBoundary(NODE3D_BDRY_TYPE);
   for (int i = 0; i < node_bdry.getSize(); i++) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(node_bdry[i].getBoundaryType() == NODE3D_BDRY_TYPE);
#endif

      int bnode_loc = node_bdry[i].getLocationIndex();

      hier::Box<3> fill_box = pgeom->getBoundaryFillBox(node_bdry[i],
                                                     interior,
                                                     gcw_to_fill);
      const hier::Index<3>& ibeg = fill_box.lower();
      const hier::Index<3>& iend = fill_box.upper();

      getcartnodebdry3d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ifirst(2), ilast(2),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         ibeg(2), iend(2),
                         ghost_cells(0), ghost_cells(1), ghost_cells(2),
                         dx,
                         bnode_loc,
                         bdry_node_conds[bnode_loc],
                         bdry_face_values.getPointer(),
                         vardata->getPointer(),
                         vardata->getDepth());

   }

}

/*
 * Function that returns the integer face boundary location 
 * corresponding to the given edge location and edge boundary
 * condition.  
 * 
 * If the edge boundary condition type or edge location are unknown, 
 * or the boundary condition type is inconsistant with the edge location 
 * an error results.
 */

int CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
   int edge_loc,
   int edge_btype)
{
   
   int ret_face = -1;
 
   switch (edge_btype) {
      case XFLOW_BC:
      case XREFLECT_BC:
      case XDIRICHLET_BC: 
      case XNEUMANN_BC: 
      {
         if (edge_loc == XLO_ZLO_3D || edge_loc == XLO_ZHI_3D ||
             edge_loc == XLO_YLO_3D || edge_loc == XLO_YHI_3D) {
            ret_face = XLO;
         } else if (edge_loc == XHI_ZLO_3D || edge_loc == XHI_ZHI_3D ||
                    edge_loc == XHI_YLO_3D || edge_loc == XHI_YHI_3D) {
            ret_face = XHI;
         }
         break; 
      }
      case YFLOW_BC:
      case YREFLECT_BC:
      case YDIRICHLET_BC: 
      case YNEUMANN_BC: 
      {
         if (edge_loc == YLO_ZLO_3D || edge_loc == YLO_ZHI_3D ||
             edge_loc == XLO_YLO_3D || edge_loc == XHI_YLO_3D) {
            ret_face = YLO;
         } else if (edge_loc == YHI_ZLO_3D || edge_loc == YHI_ZHI_3D ||
                    edge_loc == XLO_YHI_3D || edge_loc == XHI_YHI_3D) {
            ret_face = YHI;
         }
         break; 
      }
      case ZFLOW_BC:
      case ZREFLECT_BC:
      case ZDIRICHLET_BC: 
      case ZNEUMANN_BC: 
      {
         if (edge_loc == YLO_ZLO_3D || edge_loc == YHI_ZLO_3D ||
             edge_loc == XLO_ZLO_3D || edge_loc == XHI_ZLO_3D) {
            ret_face = ZLO;
         } else if (edge_loc == YLO_ZHI_3D || edge_loc == YHI_ZHI_3D ||
                    edge_loc == XLO_ZHI_3D || edge_loc == XHI_ZHI_3D) {
            ret_face = ZHI;
         }
         break; 
      }
      default: { 
         TBOX_ERROR("Unknown edge boundary condition type = "
            << edge_btype << " passed to \n" 
            << "CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry" 
            << std::endl); }
   }

   if (ret_face == -1) {
      TBOX_ERROR("Edge boundary condition type = "
         << edge_btype << " and edge location = " << edge_loc
         << "\n passed to "
         << "CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry"
         << " are inconsistant." << std::endl);
   }

   return(ret_face);

}

/*
 * Function that returns the integer face boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown, 
 * or the boundary condition type is inconsistant with the node location 
 * an error results.
 */

int CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
   int node_loc,
   int node_btype)
{

   int ret_face = -1;

   switch (node_btype) {
      case XFLOW_BC:
      case XREFLECT_BC:
      case XDIRICHLET_BC:
      case XNEUMANN_BC:
      {
         if (node_loc == XLO_YLO_ZLO || node_loc == XLO_YHI_ZLO ||
             node_loc == XLO_YLO_ZHI || node_loc == XLO_YHI_ZHI) {
            ret_face = XLO;
         } else {
            ret_face = XHI;
         }
         break;
      }
      case YFLOW_BC:
      case YREFLECT_BC:
      case YDIRICHLET_BC:
      case YNEUMANN_BC:
      {
         if (node_loc == XLO_YLO_ZLO || node_loc == XHI_YLO_ZLO ||
             node_loc == XLO_YLO_ZHI || node_loc == XHI_YLO_ZHI) {
            ret_face = YLO;
         } else {
            ret_face = YHI;
         }
         break;
      }
      case ZFLOW_BC:
      case ZREFLECT_BC:
      case ZDIRICHLET_BC:
      case ZNEUMANN_BC:
      {
         if (node_loc == XLO_YLO_ZLO || node_loc == XHI_YLO_ZLO ||
             node_loc == XLO_YHI_ZLO || node_loc == XHI_YHI_ZLO) {
            ret_face = ZLO;
         } else {
            ret_face = ZHI;
         }
        break;
      }
      default: {
         TBOX_ERROR("Unknown node boundary condition type = "
            << node_btype << " passed to \n"
            << "CartesianBoundaryUtilities3::getFaceLocationForNodeBdry"
            << std::endl); }
   }

   if (ret_face == -1) {
       TBOX_ERROR("Node boundary condition type = "
            << node_btype << " and node location = " << node_loc
            << "\n passed to "
            << "CartesianBoundaryUtilities3::getFaceLocationForNodeBdry"
            << " are inconsistant." << std::endl);
   }

   return(ret_face);

}

/*
 * Function to check 3D boundary data filling.  Arguments are:
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

int CartesianBoundaryUtilities3::checkBdryData(
   const std::string& varname,
   const hier::Patch<3>& patch,
   int data_id,
   int depth,
   const hier::IntVector<3>& gcw_to_check,
   const hier::BoundaryBox<3>& bbox,
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

   tbox::Pointer<geom::CartesianPatchGeometry<3> > pgeom =
      patch.getPatchGeometry();

   tbox::Pointer< pdat::CellData<3,double> > vardata =
      patch.getPatchData(data_id);

   std::string bdry_type_str;
   if (btype == FACE3D_BDRY_TYPE) {
      bdry_type_str = "FACE";
   } else if (btype == EDGE3D_BDRY_TYPE) {
      bdry_type_str = "EDGE";
   } else if (btype == NODE3D_BDRY_TYPE) {
      bdry_type_str = "NODE";
   } else {
      TBOX_ERROR("Unknown btype " << btype
         << " passed to CartesianBoundaryUtilities3::checkBdryData()! "
         << std::endl);
   }

   tbox::plog << "\n\nCHECKING 3D " << bdry_type_str << " BDRY DATA..." << std::endl;
   tbox::plog << "varname = " << varname << " : depth = " << depth << std::endl;
   tbox::plog << "bbox = " << bbox.getBox() << std::endl;
   tbox::plog << "btype, bloc, bcase = "
        << btype << ", = " << bloc << ", = " << bcase << std::endl;

   int idir;
   double valfact=0.0, constval=0.0, dxfact=0.0;
   int offsign;

   get3dBdryDirectionCheckValues(idir, offsign,
                                 btype, bloc, bcase);

   if (btype == FACE3D_BDRY_TYPE) {

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
            << " passed to CartesianBoundaryUtilities3::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   } else if (btype == EDGE3D_BDRY_TYPE) {

      if (bcase == XFLOW_BC || bcase == YFLOW_BC || bcase == ZFLOW_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XREFLECT_BC || bcase == YREFLECT_BC || 
                 bcase == ZREFLECT_BC) {
         valfact  = -1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XDIRICHLET_BC || bcase == YDIRICHLET_BC || 
                 bcase == ZDIRICHLET_BC) {
         valfact  = 0.0;
         constval = bstate;
         dxfact   = 0.0;
      } else if (bcase == XNEUMANN_BC || bcase == YNEUMANN_BC || 
                 bcase == ZNEUMANN_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR("Unknown bcase " << bcase
            << " passed to CartesianBoundaryUtilities3::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   } else if (btype == NODE3D_BDRY_TYPE) {

      if (bcase == XFLOW_BC || bcase == YFLOW_BC || bcase == ZFLOW_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XREFLECT_BC || bcase == YREFLECT_BC || 
                 bcase == ZREFLECT_BC) {
         valfact  = -1.0;
         constval = 0.0;
         dxfact   = 0.0;
      } else if (bcase == XDIRICHLET_BC || bcase == YDIRICHLET_BC || 
                 bcase == ZDIRICHLET_BC) {
         valfact  = 0.0;
         constval = bstate;
         dxfact   = 0.0;
      } else if (bcase == XNEUMANN_BC || bcase == YNEUMANN_BC || 
                 bcase == ZNEUMANN_BC) {
         valfact  = 1.0;
         constval = 0.0;
         dxfact   = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR("Unknown bcase " << bcase
            << " passed to CartesianBoundaryUtilities3::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   }

   hier::Box<3> gbox_to_check =
      vardata->getGhostBox() * pgeom->getBoundaryFillBox(bbox,
                                                         patch.getBox(),
                                                         gcw_to_check);

   hier::Box<3> cbox = gbox_to_check;
   hier::Box<3> dbox = gbox_to_check;
   hier::Index<3> ifirst = vardata->getBox().lower();
   hier::Index<3> ilast  = vardata->getBox().upper();

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

   pdat::CellIterator<3> id(dbox);
   for (pdat::CellIterator<3> ic(cbox); ic; ic++) {
      double checkval = valfact * (*vardata)(id(), depth) + constval;
      pdat::CellIndex<3> check = ic();
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
 * Private function to read 3D face boundary data from input database.
 */

void CartesianBoundaryUtilities3::read3dBdryFaces(
   BoundaryUtilityStrategy* bdry_strategy,
   tbox::Pointer<tbox::Database> bdry_db,
   tbox::Array<int>& face_conds,
   const hier::IntVector<3>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(bdry_strategy != (BoundaryUtilityStrategy*)NULL);
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(face_conds.getSize() == NUM_3D_FACES);
#endif

   int num_per_dirs = 0;
   for (int id = 0; id < 3; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < 3) { // face boundary input required

      for (int s = 0; s < NUM_3D_FACES; s++) {

         std::string bdry_loc_str;
         switch (s) {
            case XLO: { bdry_loc_str = "boundary_face_xlo";
                        break; }
            case XHI: { bdry_loc_str = "boundary_face_xhi";
                        break; }
            case YLO: { bdry_loc_str = "boundary_face_ylo";
                        break; }
            case YHI: { bdry_loc_str = "boundary_face_yhi";
                        break; }
            case ZLO: { bdry_loc_str = "boundary_face_zlo";
                        break; }
            case ZHI: { bdry_loc_str = "boundary_face_zhi";
                        break; }
            default: NULL_STATEMENT;
         }

         bool need_data_read = true;
         if (num_per_dirs > 0) {
            if ( periodic(0) && (s == XLO || s == XHI) ) {
               need_data_read = false;
            } else if ( periodic(1) && (s == YLO || s == YHI) ) {
               need_data_read = false;
            } else if ( periodic(2) && (s == ZLO || s == ZHI) ) {
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
                        face_conds[s] = FLOW_BC;
                     } else if (bdry_cond_str == "REFLECT") {
                        face_conds[s] = REFLECT_BC;
                     } else if (bdry_cond_str == "DIRICHLET") {
                        face_conds[s] = DIRICHLET_BC;
                        bdry_strategy->
                           readDirichletBoundaryDataEntry(bdry_loc_db,
                                                          bdry_loc_str,
                                                          s);
                     } else if (bdry_cond_str == "NEUMANN") {
                        face_conds[s] = NEUMANN_BC;
                        bdry_strategy->
                           readNeumannBoundaryDataEntry(bdry_loc_db,
                                                        bdry_loc_str,
                                                        s);
                     } else {
                     TBOX_ERROR("Unknown face boundary string = "
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

   } // if (num_per_dirs < 3)

}

/*
 * Private function to read 3D edge boundary data from input database.
 */

void CartesianBoundaryUtilities3::read3dBdryEdges(
   tbox::Pointer<tbox::Database> bdry_db,
   const tbox::Array<int>& face_conds,
   tbox::Array<int>& edge_conds,
   const hier::IntVector<3>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(face_conds.getSize() == NUM_3D_FACES);
   TBOX_ASSERT(edge_conds.getSize() == NUM_3D_EDGES);
#endif

   int num_per_dirs = 0;
   for (int id = 0; id < 3; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < 2) {  // edge boundary input required

      for (int s = 0; s < NUM_3D_EDGES; s++) {

         std::string bdry_loc_str;
         switch (s) {
            case YLO_ZLO_3D: { bdry_loc_str = "boundary_edge_ylo_zlo";
                            break; }
            case YHI_ZLO_3D: { bdry_loc_str = "boundary_edge_yhi_zlo";
                            break; }
            case YLO_ZHI_3D: { bdry_loc_str = "boundary_edge_ylo_zhi";
                            break; }
            case YHI_ZHI_3D: { bdry_loc_str = "boundary_edge_yhi_zhi";
                            break; }
            case XLO_ZLO_3D: { bdry_loc_str = "boundary_edge_xlo_zlo";
                            break; }
            case XLO_ZHI_3D: { bdry_loc_str = "boundary_edge_xlo_zhi";
                            break; }
            case XHI_ZLO_3D: { bdry_loc_str = "boundary_edge_xhi_zlo";
                            break; }
            case XHI_ZHI_3D: { bdry_loc_str = "boundary_edge_xhi_zhi";
                            break; }
            case XLO_YLO_3D: { bdry_loc_str = "boundary_edge_xlo_ylo";
                            break; }
            case XHI_YLO_3D: { bdry_loc_str = "boundary_edge_xhi_ylo";
                            break; }
            case XLO_YHI_3D: { bdry_loc_str = "boundary_edge_xlo_yhi";
                            break; }
            case XHI_YHI_3D: { bdry_loc_str = "boundary_edge_xhi_yhi";
                            break; }
            default: NULL_STATEMENT;
         }

         bool need_data_read = false;
         if (num_per_dirs == 0) {
            need_data_read = true;
         } else if ( periodic(0) && 
                     (s == YLO_ZLO_3D || s == YHI_ZLO_3D || 
                      s == YLO_ZHI_3D || s == YHI_ZHI_3D) ) { 
             need_data_read = true;
         } else if ( periodic(1) && 
                     (s == XLO_ZLO_3D || s == XLO_ZHI_3D || 
                      s == XHI_ZLO_3D || s == XHI_ZHI_3D) ) {
             need_data_read = true; 
         } else if ( periodic(2) &&
                     (s == XLO_YLO_3D || s == XHI_YLO_3D || 
                      s == XLO_YHI_3D || s == XHI_YHI_3D) ) {
             need_data_read = true;
         }

         if (need_data_read) {
            if (bdry_db->keyExists(bdry_loc_str)) {
               tbox::Pointer<tbox::Database> bdry_loc_db = 
                  bdry_db->getDatabase(bdry_loc_str);
               if (!bdry_loc_db.isNull()) {
                  if (bdry_loc_db->keyExists("boundary_condition")) {
                     std::string bdry_cond_str =
                        bdry_loc_db->getString("boundary_condition");
                     if (bdry_cond_str == "XFLOW") {
                        edge_conds[s] = XFLOW_BC;
                     } else if (bdry_cond_str == "YFLOW") {
                        edge_conds[s] = YFLOW_BC;
                     } else if (bdry_cond_str == "ZFLOW") {
                        edge_conds[s] = ZFLOW_BC;
                     } else if (bdry_cond_str == "XREFLECT") {
                        edge_conds[s] = XREFLECT_BC;
                     } else if (bdry_cond_str == "YREFLECT") {
                        edge_conds[s] = YREFLECT_BC;
                     } else if (bdry_cond_str == "ZREFLECT") {
                        edge_conds[s] = ZREFLECT_BC;
                     } else if (bdry_cond_str == "XDIRICHLET") {
                        edge_conds[s] = XDIRICHLET_BC;
                     } else if (bdry_cond_str == "YDIRICHLET") {
                        edge_conds[s] = YDIRICHLET_BC;
                     } else if (bdry_cond_str == "ZDIRICHLET") {
                        edge_conds[s] = ZDIRICHLET_BC;
                     } else if (bdry_cond_str == "XNEUMANN") {
                        edge_conds[s] = XNEUMANN_BC;
                     } else if (bdry_cond_str == "YNEUMANN") {
                        edge_conds[s] = YNEUMANN_BC;
                     } else if (bdry_cond_str == "ZNEUMANN") {
                        edge_conds[s] = ZNEUMANN_BC;
                     } else {
                        TBOX_ERROR("Unknown edge boundary string = "
                           << bdry_cond_str << " found in input." << std::endl);
                     }
   
                     bool ambiguous_type = false;
                     if (bdry_cond_str == "XFLOW" ||
                         bdry_cond_str == "XREFLECT" ||
                         bdry_cond_str == "XDIRICHLET" ||
                         bdry_cond_str == "XNEUMANN") {
                        if (s == YLO_ZLO_3D || s == YHI_ZLO_3D ||
                            s == YLO_ZHI_3D || s == YHI_ZHI_3D) {
                           ambiguous_type = true;
                        }
                     } else if (bdry_cond_str == "YFLOW" ||
                                bdry_cond_str == "YREFLECT" ||
                                bdry_cond_str == "YDIRICHLET" ||
                                bdry_cond_str == "YNEUMANN") {
                        if (s == XLO_ZLO_3D || s == XLO_ZHI_3D ||
                            s == XHI_ZLO_3D || s == XHI_ZHI_3D) {
                           ambiguous_type = true;
                        }
                     } else if (bdry_cond_str == "ZFLOW" ||
                                bdry_cond_str == "ZREFLECT" ||
                                bdry_cond_str == "ZDIRICHLET" ||
                                bdry_cond_str == "ZNEUMANN") {
                        if (s == XLO_YLO_3D || s == XHI_YLO_3D ||
                            s == XLO_YHI_3D || s == XHI_YHI_3D) {
                           ambiguous_type = true;
                        }
                     }
                     if (ambiguous_type) {
                        TBOX_ERROR("Ambiguous bdry condition " 
                           << bdry_cond_str
                           << " found for " << bdry_loc_str << std::endl);
                     }

                     std::string proper_face;
                     std::string proper_face_data;
                     bool no_face_data_found = false;
                     if (bdry_cond_str == "XFLOW" ||
                         bdry_cond_str == "XDIRICHLET" ||
                         bdry_cond_str == "XNEUMANN" ||
                         bdry_cond_str == "XREFLECT") {
                        if (s == XLO_ZLO_3D || s == XLO_ZHI_3D ||
                            s == XLO_YLO_3D || s == XLO_YHI_3D) {
                           proper_face = "XLO";
                           if (bdry_cond_str == "XFLOW" &&
                               face_conds[XLO] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "XDIRICHLET" &&
                               face_conds[XLO] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "XNEUMANN" &&
                               face_conds[XLO] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "XREFLECT" &&
                               face_conds[XLO] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "XHI";
                           if (bdry_cond_str == "XFLOW" &&
                               face_conds[XHI] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "XDIRICHLET" &&
                               face_conds[XHI] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "XNEUMANN" &&
                               face_conds[XHI] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "XREFLECT" &&
                               face_conds[XHI] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     } else if (bdry_cond_str == "YFLOW" ||
                                bdry_cond_str == "YDIRICHLET" ||
                                bdry_cond_str == "YNEUMANN" ||
                                bdry_cond_str == "YREFLECT") {
                        if (s == XLO_ZLO_3D || s == YLO_ZHI_3D ||
                            s == XLO_YLO_3D || s == XHI_YLO_3D) {
                           proper_face = "YLO";
                           if (bdry_cond_str == "YFLOW" &&
                               face_conds[YLO] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "YDIRICHLET" &&
                               face_conds[YLO] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "YNEUMANN" &&
                               face_conds[YLO] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "YREFLECT" &&
                               face_conds[YLO] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "YHI";
                           if (bdry_cond_str == "YFLOW" &&
                               face_conds[YHI] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "YDIRICHLET" &&
                               face_conds[YHI] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "YNEUMANN" &&
                               face_conds[YHI] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "YREFLECT" &&
                               face_conds[YHI] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     } else if (bdry_cond_str == "ZFLOW" ||
                                bdry_cond_str == "ZDIRICHLET" ||
                                bdry_cond_str == "ZNEUMANN" ||
                                bdry_cond_str == "ZREFLECT") {
                        if (s == XLO_ZLO_3D || s == YHI_ZLO_3D ||
                            s == XLO_ZLO_3D || s == XHI_ZLO_3D) {
                           proper_face = "ZLO";
                           if (bdry_cond_str == "ZFLOW" &&
                               face_conds[ZLO] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "ZDIRICHLET" &&
                               face_conds[ZLO] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "ZNEUMANN" &&
                               face_conds[ZLO] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "ZREFLECT" &&
                               face_conds[ZLO] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "ZHI";
                           if (bdry_cond_str == "ZFLOW" &&
                               face_conds[ZHI] != FLOW_BC) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "ZDIRICHLET" &&
                               face_conds[ZHI] != DIRICHLET_BC) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "ZNEUMANN" &&
                               face_conds[ZHI] != NEUMANN_BC) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "ZREFLECT" &&
                               face_conds[ZHI] != REFLECT_BC) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     }
                     if (no_face_data_found) {
                        TBOX_ERROR("Bdry condition " << bdry_cond_str
                           << " found for " << bdry_loc_str
                           << "\n but no " << proper_face_data 
                           << " data found for face " << proper_face << std::endl);
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
 * Private function to read 3D node boundary data from input database.
 */

void CartesianBoundaryUtilities3::read3dBdryNodes(
   tbox::Pointer<tbox::Database> bdry_db,
   const tbox::Array<int>& face_conds,
   tbox::Array<int>& node_conds,
   const hier::IntVector<3>& periodic)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!bdry_db.isNull());
   TBOX_ASSERT(face_conds.getSize() == NUM_3D_FACES);
   TBOX_ASSERT(node_conds.getSize() == NUM_3D_NODES);
#endif

   int num_per_dirs = 0;
   for (int id = 0; id < 3; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < 1) { // node boundary data required 

      for (int s = 0; s < NUM_3D_NODES; s++) {

         std::string bdry_loc_str;
         switch (s) {
            case XLO_YLO_ZLO: { bdry_loc_str = "boundary_node_xlo_ylo_zlo";
                                break; }
            case XHI_YLO_ZLO: { bdry_loc_str = "boundary_node_xhi_ylo_zlo";
                                break; }
            case XLO_YHI_ZLO: { bdry_loc_str = "boundary_node_xlo_yhi_zlo";
                                break; }
            case XHI_YHI_ZLO: { bdry_loc_str = "boundary_node_xhi_yhi_zlo";
                                break; }
            case XLO_YLO_ZHI: { bdry_loc_str = "boundary_node_xlo_ylo_zhi";
                                break; }
            case XHI_YLO_ZHI: { bdry_loc_str = "boundary_node_xhi_ylo_zhi";
                                break; }
            case XLO_YHI_ZHI: { bdry_loc_str = "boundary_node_xlo_yhi_zhi";
                                break; }
            case XHI_YHI_ZHI: { bdry_loc_str = "boundary_node_xhi_yhi_zhi";
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
                  } else if (bdry_cond_str == "ZFLOW") {
                     node_conds[s] = ZFLOW_BC;
                  } else if (bdry_cond_str == "XREFLECT") {
                     node_conds[s] = XREFLECT_BC;
                  } else if (bdry_cond_str == "YREFLECT") {
                     node_conds[s] = YREFLECT_BC;
                  } else if (bdry_cond_str == "ZREFLECT") {
                     node_conds[s] = ZREFLECT_BC;
                  } else if (bdry_cond_str == "XDIRICHLET") {
                     node_conds[s] = XDIRICHLET_BC;
                  } else if (bdry_cond_str == "YDIRICHLET") {
                     node_conds[s] = YDIRICHLET_BC;
                  } else if (bdry_cond_str == "ZDIRICHLET") {
                     node_conds[s] = ZDIRICHLET_BC;
                  } else if (bdry_cond_str == "XNEUMANN") {
                     node_conds[s] = XNEUMANN_BC;
                  } else if (bdry_cond_str == "YNEUMANN") {
                     node_conds[s] = YNEUMANN_BC;
                  } else if (bdry_cond_str == "ZNEUMANN") {
                     node_conds[s] = ZNEUMANN_BC;
                  } else {
                     TBOX_ERROR("Unknown node boundary string = "
                        << bdry_cond_str << " found in input." << std::endl);
                  }

                  std::string proper_face;
                  std::string proper_face_data;
                  bool no_face_data_found = false;
                  if (bdry_cond_str == "XFLOW" ||
                      bdry_cond_str == "XDIRICHLET" ||
                      bdry_cond_str == "XNEUMANN" ||
                      bdry_cond_str == "XREFLECT") {
                     if (s == XLO_YLO_ZLO || s == XLO_YHI_ZLO ||
                         s == XLO_YLO_ZHI || s == XLO_YHI_ZHI) {
                        proper_face = "XLO";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[XLO] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[XLO] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[XLO] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[XLO] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "XHI";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[XHI] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[XHI] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[XHI] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[XHI] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  } else if (bdry_cond_str == "YFLOW" ||
                             bdry_cond_str == "YDIRICHLET" ||
                             bdry_cond_str == "YNEUMANN" ||
                             bdry_cond_str == "YREFLECT") {
                     if (s == XLO_YLO_ZLO || s == XHI_YLO_ZLO ||
                         s == XLO_YLO_ZHI || s == XHI_YLO_ZHI) {
                        proper_face = "YLO";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[YLO] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[YLO] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[YLO] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[YLO] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "YHI";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[YHI] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[YHI] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[YHI] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[YHI] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  } else if (bdry_cond_str == "ZFLOW" ||
                             bdry_cond_str == "ZDIRICHLET" ||
                             bdry_cond_str == "ZNEUMANN" ||
                             bdry_cond_str == "ZREFLECT") {
                     if (s == XLO_YLO_ZLO || s == XHI_YLO_ZLO ||
                         s == XLO_YHI_ZLO || s == XHI_YHI_ZLO) {
                        proper_face = "ZLO";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[ZLO] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[ZLO] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[ZLO] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[ZLO] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "ZHI";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[ZHI] != FLOW_BC) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[ZHI] != DIRICHLET_BC) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[ZHI] != NEUMANN_BC) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[ZHI] != REFLECT_BC) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  }
                  if (no_face_data_found) {
                     TBOX_ERROR("Bdry condition " << bdry_cond_str
                        << " found for " << bdry_loc_str
                        << "\n but no " << proper_face_data
                        << " data found for face " << proper_face << std::endl);
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
 * 3D boundary condition checking.  Called from checkBdryData().
 */

void CartesianBoundaryUtilities3::get3dBdryDirectionCheckValues(
   int& idir, 
   int& offsign, 
   int btype, 
   int bloc, 
   int bcase)
{

   std::string bdry_type_str;

   if (btype == FACE3D_BDRY_TYPE) {

      bdry_type_str = "FACE";

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
      } else if (bloc == ZLO || bloc == ZHI) {
         idir = 2;
         if (bloc == ZLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else {
         TBOX_ERROR("Unknown boundary location " << bloc
            << " passed to CartesianBoundaryUtilities3::checkBdryData()"
            << "\n for " << bdry_type_str << " boundary " << std::endl);
      }

   } else if (btype == EDGE3D_BDRY_TYPE) {

      bdry_type_str = "EDGE";

      bool bad_case = false;
      if (bcase == XFLOW_BC || bcase == XREFLECT_BC ||
          bcase == XDIRICHLET_BC || bcase == XNEUMANN_BC) {
         idir = 0;
         if (bloc == XLO_ZLO_3D || bloc == XLO_ZHI_3D ||
             bloc == XLO_YLO_3D || bloc == XLO_YHI_3D) {
            offsign = -1;
         } else if (bloc == XHI_ZLO_3D || bloc == XHI_ZHI_3D ||
                    bloc == XHI_YLO_3D || bloc == XHI_YHI_3D) {
            offsign = 1;
         } else {
            bad_case = true;
         }
      } else if (bcase == YFLOW_BC || bcase == YREFLECT_BC ||
                 bcase == YDIRICHLET_BC || bcase == YNEUMANN_BC) {
         idir = 1;
         if (bloc == YLO_ZLO_3D || bloc == YLO_ZHI_3D ||
             bloc == XLO_YLO_3D || bloc == XHI_YLO_3D) {
            offsign = -1;
         } else if (bloc == YHI_ZLO_3D || bloc == YHI_ZHI_3D ||
                    bloc == XLO_YHI_3D || bloc == XHI_YHI_3D) {
            offsign = 1;
         } else {
            bad_case = true;
         }
      } else if (bcase == ZFLOW_BC || bcase == ZREFLECT_BC ||
                 bcase == ZDIRICHLET_BC || bcase == ZNEUMANN_BC) {
         idir = 2;
         if (bloc == YLO_ZLO_3D || bloc == YHI_ZLO_3D ||
             bloc == XLO_ZLO_3D || bloc == XHI_ZLO_3D) {
            offsign = -1;
         } else if (bloc == YLO_ZHI_3D || bloc == YHI_ZHI_3D ||
                    bloc == XLO_ZHI_3D || bloc == XHI_ZHI_3D) {
           offsign = 1;
         } else {
            bad_case = true;
         }
      }
      if (bad_case) {
         TBOX_ERROR("Unknown or ambigous bcase " << bcase
            << " passed to CartesianBoundaryUtilities3::checkBdryData()"
            << "\n for " << bdry_type_str
            << " at location " << bloc << std::endl);
      }

   } else if (btype == NODE3D_BDRY_TYPE) {

      bdry_type_str = "NODE";

      if (bcase == XFLOW_BC || bcase == XREFLECT_BC ||
          bcase == XDIRICHLET_BC || bcase == XNEUMANN_BC) {
         idir = 0;
         if (bloc == XLO_YLO_ZLO || bloc == XLO_YHI_ZLO ||
             bloc == XLO_YLO_ZHI || bloc == XLO_YHI_ZHI) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == YFLOW_BC || bcase == YREFLECT_BC ||
                 bcase == YDIRICHLET_BC || bcase == YNEUMANN_BC) {
         idir = 1;
         if (bloc == XLO_YLO_ZLO || bloc == XHI_YLO_ZLO ||
             bloc == XLO_YLO_ZHI || bloc == XHI_YLO_ZHI) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == ZFLOW_BC || bcase == ZREFLECT_BC ||
                 bcase == ZDIRICHLET_BC || bcase == ZNEUMANN_BC) {
         idir = 2;
         if (bloc == XLO_YLO_ZLO || bloc == XHI_YLO_ZLO ||
             bloc == XLO_YHI_ZLO || bloc == XHI_YHI_ZLO) {
            offsign = -1;
         } else {
           offsign = 1;
         }
      }

   } else {
      TBOX_ERROR("Unknown boundary type " << btype
         << " passed to CartesianBoundaryUtilities3::checkBdryData()"
         << "\n for " << bdry_type_str
         << " at location " << bloc << std::endl);
   }
  
}

/*
 * Private function to stuff 3D boundary contants into Fortran common blocks
 */

void CartesianBoundaryUtilities3::stuff3dBdryFortConst()
{
   stufcartbdryloc3d_(XLO, XHI, YLO, YHI, ZLO, ZHI,
                      YLO_ZLO_3D, YHI_ZLO_3D, YLO_ZHI_3D, YHI_ZHI_3D,
                      XLO_ZLO_3D, XLO_ZHI_3D, XHI_ZLO_3D, XHI_ZHI_3D,
                      XLO_YLO_3D, XHI_YLO_3D, XLO_YHI_3D, XHI_YHI_3D,
                      XLO_YLO_ZLO, XHI_YLO_ZLO, XLO_YHI_ZLO, XHI_YHI_ZLO,
                      XLO_YLO_ZHI, XHI_YLO_ZHI, XLO_YHI_ZHI, XHI_YHI_ZHI);
   stufcartbdrycond3d_(FLOW_BC,
                       XFLOW_BC, YFLOW_BC, ZFLOW_BC,
                       REFLECT_BC,
                       XREFLECT_BC, YREFLECT_BC, ZREFLECT_BC,
                       DIRICHLET_BC,
                       XDIRICHLET_BC, YDIRICHLET_BC, ZDIRICHLET_BC,
                       NEUMANN_BC,
                       XNEUMANN_BC, YNEUMANN_BC, ZNEUMANN_BC);
   s_fortran_constants_stuffed = true;
}

}
}
