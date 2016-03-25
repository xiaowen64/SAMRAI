//
// File:        EmbeddedBoundaryGeometry.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 609 $
// Modified:    $Date: 2005-09-13 15:15:49 -0700 (Tue, 13 Sep 2005) $
// Description: Compute and store geometry information about the 
//              embedded boundary
//              
// 

#ifndef included_appu_EmbeddedBoundaryGeometry_C
#define included_appu_EmbeddedBoundaryGeometry_C

#include "EmbeddedBoundaryGeometry.h"

#include "BoundaryBox.h"
#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CartesianBoundaryDefines.h"
#include "CellData.h"
#include "CellIterator.h"
#include "CellIndex.h"
#include "CellIntegerConstantRefine.h"
#include "CellDoubleConstantRefine.h"
#include "tbox/IEEE.h"
#include "IndexData.h"
#include "IntVector.h"
#include "tbox/IOStream.h"
#include "NodeData.h"
#include "NodeIterator.h"
#include "NodeIndex.h"
#include "Patch.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

// Shapes used to cut embedded boundary
#include "EmbeddedBoundaryShapePolygon.h"
#include "EmbeddedBoundaryShapeSphere.h"

#define APPU_EMBEDDED_BOUNDARY_GEOMETRY_VERSION (1)
#define EBGEOM_UNDEFINED (-1)
#define CUTCASE_UNDEFINED (-1)

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#undef USE_SINGLE_POINT_FOR_INOUT
//#define USE_SINGLE_POINT_FOR_INOUT 
//#undef USE_ARRAY_FOR_INOUT
#define USE_ARRAY_FOR_INOUT
#undef DEBUG_PRINT
//#define DEBUG_PRINT
#undef RECORD_STATS
//#define RECORD_STATS


#ifdef RECORD_STATS
#include "tbox/Statistic.h"
#include "tbox/Statistician.h"
#endif


/*
*************************************************************************
*
* External declarations for FORTRAN 77 routines
*
*************************************************************************
*/
extern "C" {
   void setebparams_(const int&, const int&, const int&, const int&,
                     const int&, const int&, const int&, const int&);

// 2D fortran subroutines
   void node2cellflag2d_(const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int&, const int&,
                         const int*,
                         int*,
                         double*);

   void setebnode2d_(const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     int*, double*);

   void setebedge2d_(const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     double*, double*);



// 3D fortran subroutines
   void node2cellflag3d_(const int&, const int&, const int&,
                         const int&, const int&, const int&,
                         const int&, const int&, const int&,
                         const int&, const int&, const int&,
                         const int*,
                         int*,
                         double*);

   void setebnode3d_(const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&, const int&,
                     const int&, const int&,
                     double*, double*);

   void setebedge3d_(const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&, const int&,
                     const int&, const int&,
                     double*, double*);

   void setebface3d_(const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&,
                     const int&, const int&, const int&,
                     const int&, const int&,
                     double*, double*);

}


   
namespace SAMRAI {
    namespace appu {

/*
*************************************************************************
*                                                                       *
*  Constructor initializes private data and reads input.                * 
*                                                                       *
*************************************************************************
*/
template<int DIM> 
EmbeddedBoundaryGeometry<DIM>::EmbeddedBoundaryGeometry(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   const tbox::Pointer<geom::CartesianGridGeometry<DIM> > grid_geom,
   const hier::IntVector<DIM>& nghosts) 
{

   d_object_name = object_name;

   /*
    * Set commons in the fortran routines used by this class.
    */
   setebparams_(SOLID,CUT,BORDER,FLOW,               // cell
                OUTSIDE,INSIDE,BOUNDARY,ONBOUNDARY); // node

   /*
    * Register with RestartManager to write restart info via
    * the "putToDatabase()" method.
    */
   tbox::RestartManager::getManager()->
      registerRestartItem(d_object_name,this);

   d_grid_geometry = grid_geom;
   
   d_cell_flag_data_id = EBGEOM_UNDEFINED;
   d_cell_vol_data_id = EBGEOM_UNDEFINED;
   d_ebdry_data_id = EBGEOM_UNDEFINED;
   d_node_flag_data_id = EBGEOM_UNDEFINED;

   /*
    * Initialize constants in the problem.
    */
   d_max_subdivides = 0;
   d_max_levels = 0;
   d_use_recursive_algs = false;
   d_compute_areas_and_normal = false;
   d_compute_cutcell_index_data = false;
   d_compute_boundary_node_data = false;

   /*
    * Initialize variables and communication algorithms.
    */
   initializeVariables(nghosts);

   /*
    * If this is a restarted case, get information from restart.
    */
   bool is_from_restart = 
      tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   /*
    * Initialize shape info.
    */
   d_use_cubes = false;
   d_use_eleven_boundary_node = false;
   d_use_eleven_inside_outside = false;
   
   /*
    * Get information from input.
    */
   if (!(input_db.isNull())) {
      getFromInput(input_db, is_from_restart);
   }

   if (d_use_cubes) {
      d_cubes_interface = new appu::CubesPatchInterface<DIM>(
         "CubesPatchInterface",
         input_db->getDatabase("CubesPatchInterface"),
         d_grid_geometry,
         nghosts);
   }

   if (d_use_eleven_boundary_node || d_use_eleven_inside_outside) {
      d_eleven_interface = new ElevenPatchInterface<DIM>(
         "ElevenPatchInterface",
         input_db->getDatabase("ElevenPatchInterface"));
   }
      
      
   /*
    * Set defaults for boundary conditions.  Initialize to bogus values
    * for error checking.
    */
   int i;
   if (DIM == 2) {
      d_edge_bdry_cond.resizeArray(NUM_2D_EDGES);
      for (i = 0; i < NUM_2D_EDGES; i++) {
         d_edge_bdry_cond[i] = tbox::IEEE::getINT_MAX();
      }
   
      d_node_bdry_cond.resizeArray(NUM_2D_NODES);
      for (i = 0; i < NUM_2D_NODES; i++) {
         d_node_bdry_cond[i] = tbox::IEEE::getINT_MAX();
      }
   }
   
   if (DIM == 3) {
      d_face_bdry_cond.resizeArray(NUM_3D_FACES);
      for (i = 0; i < NUM_3D_FACES; i++) {
         d_face_bdry_cond[i] = tbox::IEEE::getINT_MAX();
      }
      d_edge_bdry_cond.resizeArray(NUM_3D_EDGES);
      for (i = 0; i < NUM_3D_EDGES; i++) {
         d_edge_bdry_cond[i] = tbox::IEEE::getINT_MAX();
      }
      
      d_node_bdry_cond.resizeArray(NUM_3D_NODES);
      for (i = 0; i < NUM_3D_NODES; i++) {
         d_node_bdry_cond[i] = tbox::IEEE::getINT_MAX();
      }
   }

   /*
    * The boundary node stuff only works if you have at least 
    * one ghost cell.
    */
   if (d_compute_boundary_node_data) {
      for (i = 0; i < DIM; i++) {
         if (nghosts(i) <= 0) {
            TBOX_ERROR(d_object_name << ": cannot compute boundary node"
                       << " information unless at least one ghost is used."
                       << "\nAdjust the 'nghosts' argument used for "
                       << "constructing this object, or set" 
                       << "\n'compute_boundary_node_data = FALSE' in "
                       << "the input file." << endl);
         }
      }
   }


   /*
    * Timers used in this class.
    */
   t_compute_eb = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::computeEBOnLevel()");
   t_calc_node_inout = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::tagInsideOutsideNodes()");
   t_calc_volume = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::calculateVolume()");
   t_calc_area = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::calculateArea()");
   t_calc_boundary_node = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::calculateBoundaryNodes()");
   t_eb_phys_bdry = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::setEBAtPhysBoundaries()");
   t_read_write_eb = tbox::TimerManager::getManager()->
     getTimer("appu::EmbeddedBoundaryGeometry::read_writeEBto_fromFile()");
   
}
  
/*
*************************************************************************
*                                                                       *
* Empty destructor.                                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
EmbeddedBoundaryGeometry<DIM>::~EmbeddedBoundaryGeometry<DIM>()
{
}


/*
*************************************************************************
*                                                                       *
*  Initialize geometric features of embedded boundary on the supplied   *
*  patch level. This procedure generally involves three steps:          *
*                                                                       *
*    1) Set the cell flag, cell volume, and boundary cell parameters    * 
*       (e.g. volume fraction, area fractions, surface normal).  If     *
*       the supplied level is in a patch hierarchy and a coarser level  *
*       that contains embedded boundary already exists, we can use      *
*       this information to speed up the construction.                  *
*    2) Extend the embedded boundary into the physical boundary region  *
*       of the problem (e.g. set the normal appropriate to the type of  *
*       boundary condition applied).                                    *
*    3) Construct volume of cells surrounding the cut cell, which is    *
*       used to distribute conserved quantities.                        *
*                                                                       *
*  Depending on the arguments supplied, this method may be used in      *
*  one of two ways.  The first, taking only the level as an argument,   *
*  performs an exhaustive search of all cells on the level              *
*  to find the cut cells and classify the cells as FLOW, SOLID, or CUT. * 
*  The second, taking as additional arguments a hierarchy, coarser_level,
*  and possibly old_level, performs the same function but uses the      *
*  information on the coarser and old levels to narrow the search for   *
*  cut cells, making it considerably faster.  Generally, the first method 
*  is used for the coarsest level only and the second is used for all   *
*  subsequent finer levels.                                             *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::buildEmbeddedBoundaryOnLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy,
      const tbox::Pointer<hier::PatchLevel<DIM> > old_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
#endif

   /*
    * Set level ratio information required to build the 
    * embedded boundary.
    */
   setLevelRatioInformation(level);

   /*
    * Make sure the boundary, flag, and vol data have been allocated
    * on the level.
    */
   if (!(level->checkAllocated(d_ebdry_data_id))) {
      level->allocatePatchData(d_ebdry_data_id);
   }
   if (!(level->checkAllocated(d_cell_flag_data_id))) {
      level->allocatePatchData(d_cell_flag_data_id);
   }

   if (!(level->checkAllocated(d_node_flag_data_id))) {
      level->allocatePatchData(d_node_flag_data_id);
   }

   if (!(level->checkAllocated(d_cell_vol_data_id))) {
      level->allocatePatchData(d_cell_vol_data_id);
   }
   
   /*
    * Initialze flag and vol data to FLOW conditions - flag=FLOW, vol=1.0
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag = 
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,int> > node_flag = 
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > vol = 
         patch->getPatchData(d_cell_vol_data_id);

      cell_flag->fillAll(FLOW);
      node_flag->fillAll(OUTSIDE);
      vol->fillAll(1.0);
   }
   
   /*
    * If we are operating on a temporary level, or on the coarsest level
    * of the hierarchy, iterate over all cells to classify cells as
    * FLOW, SOLID, or CUT.
    *
    * If on a finer level, use the coarser level flag 
    * information to help isolate the search.  Refine the cell flag 
    * from the coarser level and operate only on cells that
    * were identified as being split cells on the coarser level.
    */

   int level_number = level->getLevelNumber();

   /*
    * Read the embedded boundary from file or generate it.  If it is read
    * from file, the grid geometry may be null since we don't need the
    * cell dx information. 
    */
   if (d_read_ebdry) {

      readLevelEmbeddedBoundaryDataFromFile(level, d_ebdry_dirname);

   } else {      

      /*
       * If the embedded boundary is to be constructed, the grid geometry
       * must be non-null and we must have a valid array of shapes from
       * input.
       */
      if (d_grid_geometry.isNull()) {
         TBOX_ERROR(d_object_name << ":buildEmbeddedBoundaryOnLevel()"
                    << "\nThe grid geometry is NULL.  Make sure the"
                    << "\ngrid geometry is supplied to the constructor."
                    << endl);
      }
         
      /*
       * Reset error parameters for cut cell computation on level.
       */
      int cut_cells_on_proc = 0;
      double l2_volume_error = 0.;
      double l2_area_error = 0.;
      double max_volume_error = 0.;
      double max_area_error = 0.;

      /*
       * Compute solid, cut, and surrounding cells on level.  
       * - If there are no coarser levels, do an exhaustive search to find
       *   the inside/outside points on the level (expensive).
       * - If there is a coarser level in the hierarchy, use the flag
       *   information from the coarser level to narrow the search for
       *   cut cells. 
       */
      t_compute_eb->start();
      if (d_use_cubes || d_use_eleven_boundary_node) {
         computeEmbeddedBoundaryOnLevelWithPackage(level,
                                                   cut_cells_on_proc);
      } else {
         computeEmbeddedBoundaryOnLevel(level,
                                        cut_cells_on_proc,
                                        l2_volume_error,
                                        l2_area_error,
                                        max_volume_error,
                                        max_area_error);
      }
      t_compute_eb->stop();
      

      /*
       * Output information about the cut cells on the level.
       */

      int number_cells = 0;
      hier::BoxArray<DIM> level_boxes = level->getBoxes();
      for (int i = 0; i < level_boxes.getNumberOfBoxes(); i++) {
         number_cells += level_boxes(i).size();
      }

      int cut_cells_on_level = tbox::MPI::sumReduction(cut_cells_on_proc);
      
      if (d_verbose) {
         
         tbox::pout << "============ Embedded Boundary on Level " << level_number 
                    << " ============"
                    << "\n  total number cells on level: " << number_cells
                    << "\n  number cut cells on level:   " << cut_cells_on_level;
         
         if (!(d_use_cubes || d_use_eleven_boundary_node)) {
            
            if (cut_cells_on_level > 0) {
               l2_volume_error = 
                  sqrt(l2_volume_error) / (double)cut_cells_on_level;
               l2_area_error = 
                  sqrt(l2_area_error) /(double)cut_cells_on_level;
            }
            tbox::pout << "\n  L2 norm of error in volume:  " 
                       << l2_volume_error;
            
            if (d_compute_areas_and_normal) {
               tbox::pout   << "\n  L2 norm of error in area:    " 
                            << l2_area_error
                            << "\n  max error in volume:         " 
                            << max_volume_error
                            << "\n  max error in area:           " 
                            << max_area_error;
            } else {
               
               tbox::pout   << "\n  AREAS NOT COMPUTED           ";
               
            }
         }
         
         
         tbox::pout << "\n===================================================="
                    << endl;
         
      }

#ifdef RECORD_STATS
      tbox::Pointer<tbox::Statistic> num_patches_l0 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberPatchesL0", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_patches_l1 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberPatchesL1", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_patches_l2 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberPatchesL2", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_patches_l3 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberPatchesL3", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_patches_l4 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberPatchesL4", "PROC_STAT");

      tbox::Pointer<tbox::Statistic> num_gridcells_l0 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberGridcellsL0", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_gridcells_l1 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberGridcellsL1", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_gridcells_l2 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberGridcellsL2", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_gridcells_l3 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberGridcellsL3", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_gridcells_l4 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberGridcellsL4", "PROC_STAT");

      tbox::Pointer<tbox::Statistic> num_cutcells_l0 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberCutcellsL0", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_cutcells_l1 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberCutcellsL1", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_cutcells_l2 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberCutcellsL2", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_cutcells_l3 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberCutcellsL3", "PROC_STAT");
      tbox::Pointer<tbox::Statistic> num_cutcells_l4 =
         tbox::Statistician::getStatistician()->
         getStatistic("NumberCutcellsL4", "PROC_STAT");

      // count number patches
      double number_patches = 0.;
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());
         number_patches += 1.0;
      }
      
      // record stats
      if (level_number == 0) {
         num_patches_l0->recordProcStat(number_patches);
         num_gridcells_l0->recordProcStat((double)number_cells);
         num_cutcells_l0->recordProcStat((double)cut_cells_on_proc);
      }
      if (level_number == 1) {
         num_patches_l1->recordProcStat(number_patches);
         num_gridcells_l1->recordProcStat((double)number_cells);
         num_cutcells_l1->recordProcStat((double)cut_cells_on_proc);
      }
      if (level_number == 2) {
         num_patches_l2->recordProcStat(number_patches);
         num_gridcells_l2->recordProcStat((double)number_cells);
         num_cutcells_l2->recordProcStat((double)cut_cells_on_proc);
      }
      if (level_number == 3) {
         num_patches_l3->recordProcStat(number_patches);
         num_gridcells_l3->recordProcStat((double)number_cells);
         num_cutcells_l3->recordProcStat((double)cut_cells_on_proc);
      }
      if (level_number == 4) {
         num_patches_l4->recordProcStat(number_patches);
         num_gridcells_l4->recordProcStat((double)number_cells);
         num_cutcells_l4->recordProcStat((double)cut_cells_on_proc);
      }
#endif

   }

   if (d_write_ebdry) {
      
      writeLevelEmbeddedBoundaryDataToFile(level, d_ebdry_dirname);

   }
   

   /*
    * Set the boundary cells on physical boundaries according to
    * the specified boundary conditions.  
    */
   setEmbeddedBoundaryAtPhysicalBoundaries(level);
   
   /*
    * Set the "surrounding volumes" for each EB cell.  This is the 
    * sum of the volumes of the cell surrounding the EB cell and is
    * used to distribute mass.
    */
   //   setSurroundingVolumes(level);   

}


/*
*************************************************************************
*                                                                       *
* Compute characteristics of the embedded boundary using Cubes.  Loop   *
* through patches and call the patch interface to Cubes.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::computeEmbeddedBoundaryOnLevelWithPackage(
   const tbox::Pointer<hier::PatchLevel<DIM> > level,
   int &cut_cells_on_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
   if (d_use_cubes) assert(!d_cubes_interface.isNull());
   if (d_use_eleven_boundary_node) assert(!d_eleven_interface.isNull());
#endif

   if (d_use_cubes && DIM != 3) {
      TBOX_ERROR(d_object_name << "::computeEmbeddedBoundaryOnLevelWithPackage"
                 << "\nCannot use Cubes - it only works for 3D cases..");
   }

#ifdef RECORD_STATS
   tbox::Pointer<tbox::Statistic> num_cutcells_patch_l0 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberCutcellsPatchL0", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_cutcells_patch_l1 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberCutcellsPatchL1", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_cutcells_patch_l2 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberCutcellsPatchL2", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_cutcells_patch_l3 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberCutcellsPatchL3", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_cutcells_patch_l4 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberCutcellsPatchL4", "PATCH_STAT");

   tbox::Pointer<tbox::Statistic> num_flowcells_patch_l0 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberFlowcellsPatchL0", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_flowcells_patch_l1 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberFlowcellsPatchL1", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_flowcells_patch_l2 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberFlowcellsPatchL2", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_flowcells_patch_l3 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberFlowcellsPatchL3", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> num_flowcells_patch_l4 =
      tbox::Statistician::getStatistician()->
      getStatistic("NumberFlowcellsPatchL4", "PATCH_STAT");

   tbox::Pointer<tbox::Timer> t_patch = tbox::TimerManager::getManager()->
      getTimer("appu::EmbeddedBoundaryGeometry::package_patch");

   tbox::Pointer<tbox::Statistic> timing_patch_l0 =
      tbox::Statistician::getStatistician()->
      getStatistic("TimingPatchL0", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> timing_patch_l1 =
      tbox::Statistician::getStatistician()->
      getStatistic("TimingPatchL1", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> timing_patch_l2 =
      tbox::Statistician::getStatistician()->
      getStatistic("TimingPatchL2", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> timing_patch_l3 =
      tbox::Statistician::getStatistician()->
      getStatistic("TimingPatchL3", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> timing_patch_l4 =
      tbox::Statistician::getStatistician()->
      getStatistic("TimingPatchL4", "PATCH_STAT");
#endif

   /*
    * Loop through patches and call the patch interface to Cubes.
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
     tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

#ifdef RECORD_STATS
     t_patch->reset();
     t_patch->start();
#endif     

     if (d_use_cubes) {
        d_cubes_interface->calculateCutCellInfo(patch,
                                                d_cell_flag_data_id,
                                                d_cell_vol_data_id,
                                                d_ebdry_data_id);
     }

     if (d_use_eleven_boundary_node) {
        d_eleven_interface->calculateCutCellInfo(patch,
                                                 d_cell_flag_data_id,
                                                 d_cell_vol_data_id,
                                                 d_node_flag_data_id,
                                                 d_ebdry_data_id);
     }


#ifdef RECORD_STATS
     t_patch->stop();

     double num_flowcells = (double)patch->getBox().size();
     double time = t_patch->getTotalWallclockTime();
     tbox::Pointer< pdat::IndexData<DIM,appu::CutCell<DIM> > > eboundary =
        patch->getPatchData(d_ebdry_data_id);
     int num_cutcells = eboundary->getNumberItems();
     num_flowcells = num_flowcells - num_cutcells;

     // record stats
     int ln = patch->getPatchLevelNumber();
     int pn = patch->getPatchNumber();
     if (ln == 0) {
        num_cutcells_patch_l0->recordPatchStat(pn,num_cutcells,d_step_count);
        num_flowcells_patch_l0->recordPatchStat(pn,num_flowcells,d_step_count);
        timing_patch_l0->recordPatchStat(pn,time,d_step_count);
     }
     if (ln == 1) {
        num_cutcells_patch_l1->recordPatchStat(pn,num_cutcells,d_step_count);
        num_flowcells_patch_l1->recordPatchStat(pn,num_flowcells,d_step_count);
        timing_patch_l1->recordPatchStat(pn,time,d_step_count);
     }
     if (ln == 2) {
        num_cutcells_patch_l2->recordPatchStat(pn,num_cutcells,d_step_count);
        num_flowcells_patch_l2->recordPatchStat(pn,num_flowcells,d_step_count);
        timing_patch_l2->recordPatchStat(pn,time,d_step_count);
     }
     if (ln == 3) {
        num_cutcells_patch_l3->recordPatchStat(pn,num_cutcells,d_step_count);
        num_flowcells_patch_l3->recordPatchStat(pn,num_flowcells,d_step_count);
        timing_patch_l3->recordPatchStat(pn,time,d_step_count);
     }
     if (ln == 4) {
        num_cutcells_patch_l4->recordPatchStat(pn,num_cutcells,d_step_count);
        num_flowcells_patch_l4->recordPatchStat(pn,num_flowcells,d_step_count);
        timing_patch_l4->recordPatchStat(pn,time,d_step_count);
     }
#endif
     
        
   }

#ifdef RECORD_STATS
   d_step_count++;
#endif

   /*
    * Count the number of cut cells
    */

   cut_cells_on_level = 0;
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());      

      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag =
         patch->getPatchData(d_cell_flag_data_id);

      const hier::Box<DIM>& interior = patch->getBox();

      for (pdat::CellIterator<DIM> ic(interior); ic; ic++) {
         pdat::CellIndex<DIM> cell(ic());
         
         if ((*cell_flag)(cell) == CUT) {
            cut_cells_on_level++;
         }
      }
   }


}

   
/*
*************************************************************************
*                                                                       *
*  Set the embedded boundary on the supplied level.  Loop over all      *
*  cells on the level and flag them as being flow, split, boundary,     *
*  or interior cells.  For split cells, we compute the volume and       *
*  area fractions.                                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::computeEmbeddedBoundaryOnLevel(
   const tbox::Pointer<hier::PatchLevel<DIM> > level,
   int &cut_cells_on_level,
   double &l2_volume_error,
   double &l2_area_error,
   double &max_volume_error,
   double &max_area_error)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
#endif

   /**************************************************************
    *  1) Set tags on the nodes then identify cut cells:
    *     - mark nodes as being "inside" or "outside" the geom
    *     - call fortran determination to set cells as "FLOW",
    *       "CUT", or "SOLID", based on the surrounding node
    *       inside/outside determination.
    **************************************************************/
   /*
    *  Set tags on nodes of the level.
    */
   t_calc_node_inout->start();
   tagInsideOutsideNodesOnLevel(level);
   t_calc_node_inout->stop();

   /*
    * Convert node tags to cell tags.
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
     tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());
     
     tbox::Pointer< pdat::NodeData<DIM,int> > node_flag =
         patch->getPatchData(d_node_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag =
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!(node_flag.isNull()));
      assert(!(cell_flag.isNull()));
      assert(!(cell_vol.isNull()));
#endif
      /*
       * Assume cell is "flow" unless computed otherwise.
       */
      cell_flag->fillAll(FLOW);
      cell_vol->fillAll(1.0);

      const hier::Box<DIM>& interior = patch->getBox();
      const hier::Index<DIM>& ifirst = interior.lower();
      const hier::Index<DIM>& ilast  = interior.upper();

      hier::IntVector<DIM> node_ghosts = node_flag->getGhostCellWidth(); 
      hier::IntVector<DIM> cell_ghosts = cell_flag->getGhostCellWidth(); 

      if (DIM == 2) {
         node2cellflag2d_(ifirst(0),ifirst(1),
                          ilast(0),ilast(1),
                          node_ghosts(0),node_ghosts(1),
                          cell_ghosts(0),cell_ghosts(1),
                          node_flag->getPointer(),
                          cell_flag->getPointer(),
                          cell_vol->getPointer());
      }
      
      if (DIM == 3) {
         node2cellflag3d_(ifirst(0),ifirst(1),ifirst(2),
                          ilast(0),ilast(1),ilast(2),
                          node_ghosts(0),node_ghosts(1),node_ghosts(2),
                          cell_ghosts(0),cell_ghosts(1),cell_ghosts(2),
                          node_flag->getPointer(),
                          cell_flag->getPointer(),
                          cell_vol->getPointer());
      }

   } // loop over patches

   /**************************************************************
    *  2) Compute volume/areas of CUT cells
    **************************************************************/
   /*
    * Set global level info - physical domain boxes, lowest index, and 
    * the lowest XYZ coordinate.
    */
   int i;
   hier::BoxArray<DIM> domain_boxes;
   domain_boxes = level->getGridGeometry()->getPhysicalDomain();
   hier::Index<DIM> domain_ilo(tbox::IEEE::getINT_MAX());
   for (int n = 0; n < domain_boxes.getNumberOfBoxes(); n++) {
      for (i = 0; i < DIM; i++) {
         if (domain_boxes(n).lower(i) < domain_ilo(i)) {
            domain_ilo(i) = domain_boxes(n).lower(i);
         }
      }
   }
   const double* domain_xlo  = d_grid_geometry->getXLower();


   /*
    * Loop over patches on level and compute volume/areas for cut cells.
    */
   cut_cells_on_level = 0;
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());      

      const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
         patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      hier::IntVector<DIM> level_ratio = pgeom->getRatio();
      
      /* 
       * Determine full volume and cell face areas.
       */
      double fullcellvol = 1.0;
      for (i = 0; i < DIM; i++) {
         fullcellvol = fullcellvol*dx[i];
      }
      
      double fullareas[2*DIM];
      for (i = 0; i < DIM; i++) {
         fullareas[2*i] = fullcellvol/dx[i];
         fullareas[2*i+1] = fullcellvol/dx[i];
      }

      double lower[DIM], upper[DIM];
      double volume_error_estimate, area_error_estimate;

      /*
       * Loop over INTERIOR cells on the patch - this can include 
       * ghost cells but not those on physical boundary - and 
       * identify cut cells.  Embedded boundary cells that lie on 
       * the physical boundary require special treatment, performed
       * in the "setEmbeddedBoundaryAtPhysicalBoundaries()" method.
       */

      const hier::Box<DIM>& interior = patch->getBox();

      tbox::Pointer< pdat::IndexData<DIM,appu::CutCell<DIM> > > eboundary =
         patch->getPatchData(d_ebdry_data_id);
      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag =
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);
      tbox::Pointer< pdat::NodeData<DIM,int> > node_flag =
         patch->getPatchData(d_node_flag_data_id);
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!(eboundary.isNull()));
      assert(!(cell_flag.isNull()));
      assert(!(cell_vol.isNull()));
      assert(!(node_flag.isNull()));
#endif         

      for (pdat::CellIterator<DIM> ic(interior); ic; ic++) {
         pdat::CellIndex<DIM> cell(ic());
         
         if ((*cell_flag)(cell) == CUT) {

           /*
            * Compute the volume and areas.  
            *
            * Fast approximation:
            *   If the cell is not going to be sub-divided, the areas
            *   and normal are not needed (d_compute_areas_and_normal = false),
            *   and we do not need to keep a list of cut cells, do a quick
            *   estimation of the volume and exit. 
            *
            * Standard computation:
            *   - create a CutCell data struct
            *   - compute the information for the struct
            *   - append the struct to the list of IndexData<DIM> for the patch
            */
           if (d_max_subdivides == 0 && 
               d_compute_areas_and_normal &&
               !d_compute_cutcell_index_data) {

              (*cell_vol)(cell) = 0.5; 
              volume_error_estimate = 0.25 * fullcellvol;

           } else {

              /*
               * Create a cut cell.
               */
              appu::CutCell<DIM> cc(ic());
                  
              /*
               * Compute upper and lower extents of the cell
               */
              for (i = 0; i < DIM; i++) {
                 lower[i]= domain_xlo[i] + (cell(i)-domain_ilo(i))*dx[i];
                 upper[i]= domain_xlo[i] + (cell(i)-domain_ilo(i)+1)*dx[i];
              }         
                  
              /*
               * Compute volume, areas, normal, etc. for the cut cell
               */
              volume_error_estimate = 0.;
              area_error_estimate = 0.;
              calculateCutCellInformation(cc,
                                          lower,
                                          upper,
                                          fullcellvol,
                                          fullareas,
                                          volume_error_estimate,
                                          area_error_estimate);
               
              /*
               * Compute boundary node information, if desired. 
               */
              if (d_compute_boundary_node_data) {

                 t_calc_boundary_node->start();
                 calculateBoundaryNodeInformation(cc,
                                                  node_flag,
                                                  lower,
                                                  upper,
                                                  dx);
                 t_calc_boundary_node->stop();

              }
                    
                        
              /*
               * Add the cut cell to the list of index data for the patch
               */
              eboundary->appendItem(ic(), cc);
              (*cell_vol)(cell) = cc.getVolume();                  
                  
           }

            /*
             * Accumulate information on errors in volume fraction
             * computation on the cut cells.
             */
            cut_cells_on_level++;
            l2_volume_error += 
               volume_error_estimate * volume_error_estimate;
            l2_area_error += 
               area_error_estimate * area_error_estimate;
            max_volume_error = 
               tbox::Utilities::dmax(volume_error_estimate,
                                     max_volume_error);
            max_area_error = 
               tbox::Utilities::dmax(area_error_estimate,
                                     max_area_error);


                      
                  

         } // if cell is cut
         
      } // cell iterator
      
#ifdef DEBUG_PRINT
      tbox::plog << "----computeEmbeddedBoundaryOnLevel() "  << endl;
      tbox::plog << "Patch: " << patch->getPatchNumber() 
                 << "  Box: " << patch->getBox() << endl;
      tbox::plog << "Cell Flag:" << endl;
      cell_flag->print(cell_flag->getGhostBox());
      tbox::plog << "Node Flag:" << endl;
      node_flag->print(node_flag->getGhostBox());
#endif

   } // loop over patches
   
   
}


/*
*************************************************************************
*                                                                       *
* Classify all nodes on the level as being INSIDE or OUTSIDE the        *
* geometry.                                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::tagInsideOutsideNodesOnLevel(
   const tbox::Pointer<hier::PatchLevel<DIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
#endif

   int i;

   /*
    * Loop over patches on level.
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());
      
      tbox::Pointer< pdat::NodeData<DIM,int> > node_flag =
         patch->getPatchData(d_node_flag_data_id);
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!(node_flag.isNull()));
#endif
      
      const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
         patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlower = pgeom->getXLower();
      hier::IntVector<DIM> level_ratio = pgeom->getRatio();
      
      const hier::Box<DIM> patch_box = patch->getBox();
      
      /*
       * Assume node is "outside" unless computed otherwise.
       */
      node_flag->fillAll(OUTSIDE);

#ifdef USE_ARRAY_FOR_INOUT
      int nx[DIM];
      double origin[DIM];
      // NOTE: we are finding NODE locations.
      for (i = 0; i < DIM; i++) {
         nx[i] = patch->getBox().numberCells(i) + 2*d_ebdry_nghosts(i) + 1;
         origin[i] = xlower[i] - (double)d_ebdry_nghosts(i)*dx[i];
      }
      
      if (d_use_eleven_inside_outside) {
         d_eleven_interface->isInside(nx,
                                      dx,
                                      origin,
                                      node_flag->getPointer());
      } else {
         
         /*
          * Native SAMRAI-defined shapes
          */
         doNativeShapeInsideOutside(nx,dx,origin,node_flag->getPointer());
                        
      }
#endif

#ifdef USE_SINGLE_POINT_FOR_INOUT
      double node_xyz[DIM];

      /*
       * Loop over nodes on the patch.
       */
      for (typename pdat::NodeIterator<DIM> ni(patch_box); ni; ni++) {

         /*
          * Compute upper and lower extents of the cell and classify
          * it.
          */
         pdat::NodeIndex<DIM> node = ni();
         for (i = 0; i < DIM; i++) {
           node_xyz[i]= xlower[i] + (double)node(i)*dx[i];
         }

         /*
          * determine if node is "inside"
          */
         bool is_inside = false;         
         if (d_use_eleven_inside_outside) {

            is_inside = d_eleven_interface->isInside(node_xyz);

         } else {
            
            /*
             * Native SAMRAI-defined shapes
             */
            
            for (n = 0; n < d_shapes.getSize(); n++) {
               is_inside = d_shapes[n]->isInside(node_xyz);
               if (is_inside) {
                  (*node_flag)(node) = INSIDE;
                  break;
               }
            }
         }
         
      } // iteration over nodes on patch
      
#endif

   } // loop over patches of level
   
}

/*
*************************************************************************
*	                                                                *
* Compute the total volume on the supplied level.                       *
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
EmbeddedBoundaryGeometry<DIM>::computeTotalVolumeOnLevel(
   const tbox::Pointer<hier::PatchLevel<DIM> > level)
{
   double vol_this_proc = 0.;

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());
      
      const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
         patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      
      /*
       * Compute actual cell volume.
       */
      double fullcellvol = 1.0;
      for (int i = 0; i < DIM; i++) {
         fullcellvol = fullcellvol*dx[i];
      }
   
      
      /*
       * Access the volume fraction.
       */
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);

      /*
       * Multiply volume fraction times the actual cell volume and 
       * add to sum to compute the volume amongst patches located on 
       * this processor.
       */
      for (pdat::CellIterator<DIM> ic(patch->getBox()); ic; ic++) {
         vol_this_proc += (*cell_vol)(ic()) * fullcellvol;
      }

   } // loop over patches

   /*
    * Compute the total volume by summing processor contributions.
    */
   double total_volume = tbox::MPI::sumReduction(vol_this_proc);

   return(total_volume);
}



/*
*************************************************************************
*                                                                       *
* Register VisIt data writer to write data to plot files that may       *
* be postprocessed by the VisIt tool.                                   *
*                                                                       *
*************************************************************************
*/
#ifdef HAVE_HDF5
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::registerVisItDataWriter(
   tbox::Pointer<appu::VisItDataWriter<DIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(visit_writer.isNull()));
#endif

   d_visit_writer = visit_writer;

   /*
    * Register the object that supplies functions to write out materials
    * data.  In this case, "this" object does it.
    */
   d_visit_writer->setMaterialsDataWriter(this);
   
   /*
    * If the cell flag and volume fraction variables have been created
    * register them with the visit data writer.
    */
   if (d_cell_flag_data_id > 0) {
      d_visit_writer->registerPlotQuantity("Cell Flag",
                                           "SCALAR",
                                           d_cell_flag_data_id);
   }

   if (d_cell_vol_data_id > 0) {
      d_visit_writer->registerPlotQuantity("Volume Frac",
                                           "SCALAR",
                                           d_cell_vol_data_id);
   }

   if (d_node_flag_data_id > 0) {
      d_visit_writer->registerPlotQuantity("Node InOut Flag",
                                           "SCALAR",
                                           d_node_flag_data_id);
   }


   /*
    * Register flow and solid as a "material" data types.
    */
   tbox::Array<string> names(2);
   names[0] = "Flow";
   names[1] = "Solid";
   d_visit_writer->registerMaterialNames(names);
   
   
}
#endif

/*
*************************************************************************
*	                                                                *
* Return the descriptor ID for the patch data which holds the following *
*                                                                       *
*   - cell flag data id - integer at each cell, specifying it as solid, *
*     cut, border, or flow.                                             *
*                                                                       *
*   - cell volume id - double at each cell specifying volume fraction.  *
*                                                                       *
*   - cut cell data id - IndexData<DIM,appu::CutCell<DIM>> at selected cells     *
*     holding information like the volume, area fraction,               *
*     normal, etc.                                                      *
*                                                                       *
*   - node flag id - int at each node specifying whether node is       *
*     flagged as inside or outside the embedded boundary                *
*                                                                       *
*************************************************************************
*/
template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::getCellFlagDataId() const
{
   return(d_cell_flag_data_id);
}

template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::getCellVolumeDataId() const
{
   return(d_cell_vol_data_id);
}

template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::getIndexCutCellDataId() const
{
   return(d_ebdry_data_id);
}

template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::getNodeInsideOutsideDataId() const
{
   return(d_node_flag_data_id);
}


/*
*************************************************************************
*                                                                       *
* Pack "volume" data (Vis materials plot quantity) into the supplied    *
* double precision buffer.                                              *
*                                                                       *
* The returned integer value allows savings of storage by returning     *
* whether the volume of all cells on the patch is 1.0, or all 0.0,      *
* or mixed.                                                             *
*  appu::VisMaterialsDataStrategy<DIM>::VISIT_ALL_ONE  - all cells volume 1. *
*  appu::VisMaterialsDataStrategy<DIM>::VISIT_ALL_ZERO - all cells volume 0. *
*  appu::VisMaterialsDataStrategy<DIM>::VISIT_MIXED    - data is packed      *
*                                                                       *
*************************************************************************
*/

template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::packMaterialFractionsIntoDoubleBuffer(
   double* dbuffer,
   const hier::Patch<DIM>& patch,
   const hier::Box<DIM>& region,
   const string& material_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((region * patch.getBox()) == region);
#endif

   tbox::Pointer< pdat::CellData<DIM,double> > cell_vol  =
      patch.getPatchData(d_cell_vol_data_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!cell_vol.isNull());
#endif

   const hier::Box<DIM>& data_box = cell_vol->getGhostBox();
   const int box_w0 = region.numberCells(0);
   const int dat_w0 = data_box.numberCells(0);
   const int box_w1 = region.numberCells(1);
   int dat_w1 = 0;
   int box_w2 = 1;
   if (DIM == 3) {
      dat_w1 = data_box.numberCells(1);
      box_w2 = region.numberCells(2);
   }

   if (material_name == "Flow") {
      const double *const vol = cell_vol->getPointer();
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

      for (int i2 = 0; i2 < box_w2; i2++) {
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               dbuffer[buf_b1+i0] = vol[dat_indx];
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
         dat_b2 += dat_w1 * dat_w0;
      }
   }


   if (material_name == "Solid") {
      const double *const vol = cell_vol->getPointer();
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());
      double solid_frac = 0.;

      for (int i2 = 0; i2 < box_w2; i2++) {
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               solid_frac = 1.0 - vol[dat_indx];
               dbuffer[buf_b1+i0] = tbox::Utilities::dabs(solid_frac);       
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
         dat_b2 += dat_w1 * dat_w0;
      }
   }

   int return_val = appu::VisMaterialsDataStrategy<DIM>::VISIT_MIXED;
   return(return_val);
}





      
/*
*************************************************************************
*                                                                       *
* Apply the operations to construct embedded boundary on fine level     *
* wherever the cell tag indicates it being a cut cell (this information *
* is provided by the coarser level.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::postprocessRefine(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const hier::Box<DIM>& fine_box,
   const hier::IntVector<DIM>& ratio)
{
   tbox::Pointer< pdat::IndexData<DIM,appu::CutCell<DIM> > > feboundary =
      fine.getPatchData(d_ebdry_data_id);
   tbox::Pointer< pdat::CellData<DIM,int> > fcell_flag =
      fine.getPatchData(d_cell_flag_data_id);
   tbox::Pointer< pdat::CellData<DIM,double> > fcell_vol =
      fine.getPatchData(d_cell_vol_data_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(fcell_flag.isNull()));
   assert(!(fcell_vol.isNull()));
#endif

   const hier::Box<DIM> coarse_box = hier::Box<DIM>::coarsen(fine_box, ratio);
   
   hier::Box<DIM> fine_patch_box = fine.getBox();
   const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
      fine.getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* domain_xlo  = d_grid_geometry->getXLower();
   hier::IntVector<DIM> level_ratio = pgeom->getRatio();
   
   /* 
    * Determine full volume and cell face areas for fine cells.
    */
   int i;
   double fullcellvol = 1.0;
   for (i = 0; i < DIM; i++) {
      fullcellvol = fullcellvol*dx[i];
   }
   
   double fullareas[2*DIM];
   for (i = 0; i < DIM; i++) {
      fullareas[2*i] = fullcellvol/dx[i];
      fullareas[2*i+1] = fullcellvol/dx[i];
   }
   
   int cell_type;
   double lower[DIM], upper[DIM];
   double volume_error_estimate, area_error_estimate;
   
   /*
    * Loop over INTERIOR cells on the patch - this can include 
    * ghost cells but not those on physical boundary - and 
    * identify cut cells.  Embedded boundary cells that lie on 
    * the physical boundary require special treatment, performed
    * in the "setEmbeddedBoundaryAtPhysicalBoundaries()" method.
    */
   pdat::CellIndex<DIM> fine_cell;
   pdat::CellIndex<DIM> coarse_cell;
   
   for (pdat::CellIterator<DIM> ic(fine_box); ic; ic++) {
      
      fine_cell = ic();
      int fine_cell_flag = (*fcell_flag)(fine_cell);
      
      if (fine_cell_flag == BORDER) {
         
         /*
          * The refine operation interpolated "Border" cells on the 
          * coarse level to the fine level. These border cells will 
          * ALWAYS be "flow" cells on the fine level, so tag them as 
          * such here.
          */
         (*fcell_flag)(fine_cell) = FLOW;
         (*fcell_vol)(fine_cell) = 1.0;
         
      } else if (fine_cell_flag == CUT) {

         /*
          * "Cut" cells on the coarse level may be solid, cut, border,
          * or flow cells on the fine. Go through each fine cell and 
          * determine which.
          */
         coarse_cell = fine_cell/ratio;
         
         /*
          * Initialize cell as being "flow" and check if otherwise.
          */
         (*fcell_flag)(fine_cell) = FLOW;
         (*fcell_vol)(fine_cell) = 1.0;
         
         /*
          * Compute upper and lower extents of the cell and classify
          * it.
          */
         for (i=0; i < DIM; i++) {
            lower[i]= domain_xlo[i] + fine_cell(i)*dx[i];
            upper[i]= domain_xlo[i] + (fine_cell(i)+1)*dx[i];
         }
         cell_type = classifyCell(lower, upper);
         
         if (cell_type == CUT) {
            appu::CutCell<DIM> bc(fine_cell);

            volume_error_estimate = 0.;
            area_error_estimate = 0.;
            calculateCutCellInformation(bc,
                                        lower,
                                        upper,
                                        fullcellvol,
                                        fullareas,
                                        volume_error_estimate,
                                        area_error_estimate);
  
            feboundary->appendItem(ic(), bc);
            
            (*fcell_flag)(fine_cell) = CUT;
            (*fcell_vol)(fine_cell) = bc.getVolume();

            /*
             * Accumulate information on errors in volume fraction
             * computation on the cut cells.
             * NOTE:  If we decide to go back and use the refine operations
             * we'll have to re-institute this
             */
#if 0
            cut_cells_on_level++;
            l2_volume_error += 
               volume_error_estimate * volume_error_estimate;
            l2_area_error += 
               area_error_estimate * area_error_estimate;
            max_volume_error = 
               tbox::Utilities::dmax(volume_error_estimate,
                                     max_volume_error);
            max_area_error = 
               tbox::Utilities::dmax(area_error_estimate,
                                     max_area_error);
#endif
            
               
         } else if (cell_type == SOLID) {
            
            (*fcell_flag)(fine_cell) = SOLID;
            (*fcell_vol)(fine_cell) = 0.0;
            
         }
         
      } // cell was cut cell (on coarse level)

   } // iterate over cells of fine box
   
}

/*
*************************************************************************
*	                                                                *
* Unused (needed because we supply postprocessRefine())                 *
*                                                                       *
*************************************************************************
*/
template<int DIM> hier::IntVector<DIM> 
EmbeddedBoundaryGeometry<DIM>::getRefineOpStencilWidth() const
{
   return(hier::IntVector<DIM>(0));
}

/*
*************************************************************************
*	                                                                *
* Unused (needed because we supply postprocessRefine())                 *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::setPhysicalBoundaryConditions(
   hier::Patch<DIM>& patch,
   const double fill_time,
   const hier::IntVector<DIM>& ghost_width_to_fill)
{
   (void) patch;
   (void) fill_time;
   (void) ghost_width_to_fill;
}

/*
*************************************************************************
*                                                                       *
* Write data to restart.                                                *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::putToDatabase( 
   tbox::Pointer<tbox::Database> db)
{
   db->putInteger("APPU_EMBEDDED_BOUNDARY_GEOMETRY_VERSION",
                 APPU_EMBEDDED_BOUNDARY_GEOMETRY_VERSION);
 
   db->putString("d_flow_type", d_flow_type);

   tbox::Array<int> tmp_array(DIM);
   for (int i = 0; i < DIM; i++) {
      tmp_array[i] = d_ebdry_nghosts(i);
   }
   db->putIntegerArray("d_ebdry_nghosts", tmp_array);
}


/*
*************************************************************************
*                                                                       *
* Construct Variables and communication schedules used to define the    *
* embedded boundary and its inter-level operations.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::initializeVariables(
   const hier::IntVector<DIM>& nghosts)
{
   d_ebdry_nghosts = nghosts;
   
   /*
    * Initialize variables used to define the embedded boundary
    */
   hier::VariableDatabase<DIM>* variable_db = hier::VariableDatabase<DIM>::getDatabase();

   d_ebdry_var = 
     new pdat::IndexVariable<DIM,appu::CutCell<DIM> >("Embedded Boundary");
   d_cell_flag_var = 
     new pdat::CellVariable<DIM,int>("Cell Flag",1);
   d_cell_vol_var = 
     new pdat::CellVariable<DIM,double>("Cell Volume Fraction",1);
   d_node_flag_var = 
     new pdat::NodeVariable<DIM,int>("Node Flag",1);
   
   tbox::Pointer<hier::VariableContext> current =
      variable_db->getContext("CURRENT");
   tbox::Pointer<hier::VariableContext> scratch =
      variable_db->getContext("SCRATCH");
   
   d_ebdry_data_id = 
     variable_db->registerVariableAndContext(d_ebdry_var,
                                             current,
                                             nghosts);

   d_cell_flag_data_id = 
     variable_db->registerVariableAndContext(d_cell_flag_var,
                                             current,
                                             nghosts);

   d_cell_vol_data_id = 
     variable_db->registerVariableAndContext(d_cell_vol_var,
                                             current,
                                             nghosts);

   d_node_flag_data_id = 
     variable_db->registerVariableAndContext(d_node_flag_var,
                                             current,
                                             nghosts);

   d_ebdry_scratch_id = 
     variable_db->registerVariableAndContext(d_ebdry_var,
                                             current,
                                             nghosts);

   d_cell_flag_scratch_id = 
     variable_db->registerVariableAndContext(d_cell_flag_var,
                                             current,
                                             nghosts);

   d_cell_vol_scratch_id = 
     variable_db->registerVariableAndContext(d_cell_vol_var,
                                             current,
                                             nghosts);
       

   /*
    * Create a refine algorithm to refine the flag data, volume
    * fraction data, and embedded boundary.  
    */
   d_ebdry_refine_alg = new xfer::RefineAlgorithm<DIM>();

   tbox::Pointer<pdat::CellIntegerConstantRefine<DIM> > const_refine_int =
      new pdat::CellIntegerConstantRefine<DIM>();

   
   /*
    * The flag data and volume fraction data both use CONSTANT_REFINE.  
    * Since we store the volume *fraction*, not the cell volume, 
    * constant refine is the appropriate refine operator.
    */

   d_ebdry_refine_alg->registerRefine(d_cell_flag_data_id, //dst
                                      d_cell_flag_data_id, //src
                                      d_cell_flag_scratch_id, //scratch
                                      const_refine_int);

   tbox::Pointer<pdat::CellDoubleConstantRefine<DIM> > const_refine_double =
      new pdat::CellDoubleConstantRefine<DIM>();
   

   d_ebdry_refine_alg->registerRefine(d_cell_vol_data_id, //dst
                                      d_cell_vol_data_id, //src
                                      d_cell_vol_scratch_id, //scratch
                                      const_refine_double);

   /*
    * The embedded boundary uses the postprocessRefine() method, but
    * does not use any other refine operations.  Hence, the refine op is
    * set to NULL.
    */
   d_ebdry_refine_alg->registerRefine(d_ebdry_data_id, //dst
                                      d_ebdry_data_id, //src
                                      d_ebdry_data_id, //scratch
                                      NULL);

       
}


/*
*************************************************************************
*                                                                       *
* Private function that makes the necessary calls to compute volume,    *
* area fraction, normals, etc. on cut cells.                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::calculateCutCellInformation(
   appu::CutCell<DIM>& cut_cell,
   const double *lower,
   const double *upper,
   const double& fullcellvol,
   const double *fullareas,
   double& volume_error_estimate,
   double& area_error_estimate)
{

   double volume;
   int subdivide_level = 0;
   volume_error_estimate = 0.;

   /*
    * We can calculate the cut-cell volume and area information using
    * recursive or non-recursive algorithms.  See the header descriptions
    * for a discussion of each approach.
    */
   t_calc_volume->start();
   if (d_use_recursive_algs) {
      
      int return_val = recursiveCalculateVolume(lower,
                                                upper,
                                                volume,
                                                volume_error_estimate,
                                                subdivide_level);
      // approximate volume,error if we reach max subdivisions
      if (return_val == 1) {
         volume = 0.5 * fullcellvol;
         volume_error_estimate = 0.25 * fullcellvol;
      }

   } else {  // non-recursive
      
      volume = calculateVolume(lower,
                               upper,
                               volume_error_estimate);
   }
   t_calc_volume->stop();
      
   cut_cell.setVolume(volume/fullcellvol);  
   volume_error_estimate = volume_error_estimate/fullcellvol;

   if (d_compute_areas_and_normal) {
      
      t_calc_area->start();

      double areas[2*DIM];
      area_error_estimate = 1.;

      for (int f = 0; f < 2*DIM; f++) {
         subdivide_level = 0;
         double face_area_error_estimate = 0.;

         if (d_use_recursive_algs) {

            int return_val = recursiveCalculateArea(lower,
                                                    upper,
                                                    f,
                                                    areas[f],
                                                    face_area_error_estimate,
                                                    subdivide_level);
            // approximate area,error if we reach max subdivisions
            if (return_val == 1) {
               areas[f] = 0.5 * fullareas[f];
               face_area_error_estimate = 0.25 * fullareas[f];
            }

         } else { // non-recursive

            areas[f] = calculateArea(lower,
                                     upper,
                                     f,
                                     face_area_error_estimate);
         }
         
         cut_cell.setArea(areas[f]/fullareas[f],f);
         face_area_error_estimate = face_area_error_estimate/fullareas[f];
         area_error_estimate += face_area_error_estimate/((double)2*DIM);
      }
               
      t_calc_area->stop();

      int i;
      double areafront = 0.0;
      double areavector[DIM];
      for (i = 0; i < DIM; i++) {
         areavector[i] = areas[2*i] - areas[2*i+1];
         areafront = areafront + areavector[i]*areavector[i];
      }
      areafront = sqrt(areafront);
      cut_cell.setFrontArea(areafront);

      /*
       * ANDY, FIX THIS
       */
      for (i = 0; i < DIM; i++) {
         cut_cell.setFrontCentroid(0.,i);
      }
      

      /*
       * Avoid problems with normal computation if frontal area is very small
       * by resetting areafront to be large.  This will render the normal to
       * be zero.
       */
      
      if (tbox::Utilities::deq(areafront, 0.0)) {
         areafront = tbox::IEEE::getDBL_MAX();
      }
      
      double tnormal[DIM];
      for (i =0; i < DIM; i++) {
         tnormal[i] = areavector[i]/areafront;
         cut_cell.setNormal(tnormal[i], i);
      }
      cut_cell.setNewBase(tnormal);
   }

}

/*
*************************************************************************
*                                                                       *
* Private function that computes the boundary nodes and frontal 
* centroid on the cut cell.                                             *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::calculateBoundaryNodeInformation(
   appu::CutCell<DIM>& cut_cell,
   tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
   const double *lower,
   const double *upper,
   const double *dx)
{

   if (!(CutCell<DIM>::boundaryNodesEnabled())) {
      TBOX_ERROR(d_object_name << ":calculateBoundaryNodeInformation()"
                 << "\nBoundary node data cannot be computed"
                 << "\nbecause storage on CutCell is not enabled."
                 << "\nEither set 'compute_boundary_node_data' to"
                 << "\nfalse in input or use the static function"
                 << "\nappu::CutCell<DIM>::enableBoundaryNodeStorage()"
                 << "\nto enable storage." << endl);
   }
   
   hier::Index<DIM> ic = cut_cell.getIndex();

   int i, n, ilo, ihi;

   /* 
    * Loop over "edges" in each direction to mark the boundary nodes.
    */
   // set nedges = 2^(DIM-1)
   int nedges = 1;  
   for (i = 0; i < DIM-1; i++) {
      nedges *= 2;
   }

   for (i = 0; i < DIM; i++) {
      for (n = 0; n < nedges; n++) {
         hier::IntVector<DIM> nlo(0);
         hier::IntVector<DIM> nhi(0);
         if ((i == 0) && (n == 0)) {
            // nlo = 0,0, nhi = 1,0 
            nhi(0) = 1;
            ilo = 0;
            ihi = 1;
         } else if ((i == 0) && (n == 1)) {
            // nlo = 0,1, nhi = 1,1
            nlo(1) = 1;
            nhi(0) = 1;
            nhi(1) = 1;
            ilo = 2;
            ihi = 3;
         }
         if ((i == 1) && (n == 0)) {
            // nlo = 0,0  nhi = 0,1
            nhi(1) = 1;
            ilo = 0;
            ihi = 2;
         } else if ((i == 1) && (n == 1)) {
            // nlo = 1,0  nhi = 1,1
            nlo(0) = 1;
            nhi(0) = 1;
            nhi(1) = 1;
            ilo = 1;
            ihi = 3;
         }
         pdat::NodeIndex<DIM> nodelo(ic,nlo);
         pdat::NodeIndex<DIM> nodehi(ic,nhi);
         if (((*node_flag)(nodelo) == INSIDE) && 
             ((*node_flag)(nodehi) == OUTSIDE)) {
            (*node_flag)(nodelo) = BOUNDARY;
         }
         if (((*node_flag)(nodelo) == OUTSIDE) && 
             ((*node_flag)(nodehi) == INSIDE)) {
            (*node_flag)(nodehi) = BOUNDARY;
         } 
         
      } // loop over n
   } // loop over i

   /*
    * Set the boundary node information for the cut cell
    */
   hier::Box<DIM> one_cell_box(ic,ic);
   for (typename pdat::NodeIterator<DIM> ni(one_cell_box); ni; ni++) {
      pdat::NodeIndex<DIM> boundary_node = ni();
      if ((*node_flag)(boundary_node) == BOUNDARY || 
          (*node_flag)(boundary_node) == ONBOUNDARY) {

         BoundaryNode<DIM> bn(boundary_node);
         if ((*node_flag)(boundary_node) == ONBOUNDARY) {
            bn.setNodeOnBoundary();
         }

         bn.setNearestNeighborNodes(node_flag,ic);
         
         cut_cell.setBoundaryNode(bn);

      } // if node_type == BOUNDARY
   } // iterator over nodes of the cut cell

   /*
    * Set the frontal centroid.  This only works if areas are computed.
    */
   if (d_compute_areas_and_normal) {
      
      double midpt[DIM], centroid[DIM], node_loc[DIM];
      int cut_case = CUTCASE_UNDEFINED;
      int flow_side[2]; // 0 = lower, 1 = upper
      const double* area_frac = cut_cell.getArea();

      double fullcellvol = 1.0;
      for (i = 0; i < DIM; i++) {
         fullcellvol = fullcellvol*dx[i];
      }
      double fullareas[2*DIM], areas[2*DIM];
      for (i = 0; i < DIM; i++) {
         fullareas[2*i] = fullcellvol/dx[i];
         fullareas[2*i+1] = fullcellvol/dx[i];
         areas[2*i] = area_frac[2*i] * fullareas[2*i];
         areas[2*i+1] = area_frac[2*i+1] * fullareas[2*i+1];
         
      }

      if (DIM == 2) {      

         /*
          * Identify cut_case:  
          * 0 - both X faces cut
          * 1 - both Y faces cut
          * 2 - Xlo/Ylo cut
          * 3 - Xlo/Yhi cut
          * 4 - Xhi/Ylo cut
          * 5 - Xhi/Yhi cut
          */
         if ((areas[0] < dx[DIM-1]) && (areas[DIM-1] < dx[DIM-1])) {
            cut_case = 0;
         } else if ((areas[DIM] < dx[0]) && (areas[2*DIM-1] < dx[0])) {
            cut_case = 1;
         } else if ((areas[0] < dx[DIM-1]) && (areas[DIM] < dx[0])) {
            cut_case = 2;
         } else if ((areas[0] < dx[DIM-1]) && (areas[2*DIM-1] < dx[0])) {
            cut_case = 3;
         } else if ((areas[1] < dx[DIM-1]) && (areas[DIM] < dx[0])) {
            cut_case = 4;
         } else if ((areas[1] < dx[DIM-1]) && (areas[2*DIM-1] < dx[0])) {
            cut_case = 5;
         } else {
            cut_case = CUTCASE_UNDEFINED;
         }         

         /*
          * The CUTCASE_UNDEFINED can occur if one of the nodes of the 
          * cut cell lies on the boundary, in which case all the areas 
          * may be 1.0 or 0.0 so it is impossible to identify the case.  
          * For these situations, assume the centroid should not be used 
          * for any calculations so set it to a NaN value.
          */
         
         if (cut_case == 0) {
            flow_side[0] = 0;
            if (tbox::Utilities::deq(areas[DIM],0.)) {
               flow_side[1] = 0;
            } else {
               flow_side[1] = 1;
            }
         }
         if (cut_case == 1) {
            if (tbox::Utilities::deq(areas[0],0.)) {
               flow_side[0] = 0;
            } else {
               flow_side[0] = 1;
            }
            flow_side[1] = 0;
         }
         if (cut_case == 2) {
            if (tbox::Utilities::deq(areas[1],0.)) {
               flow_side[0] = 0;
            } else {
               flow_side[0] = 1;
            }
            if (tbox::Utilities::deq(areas[2*DIM-1],0.)) {
               flow_side[1] = 0;
            } else {
               flow_side[1] = 1;
            }
         }
         if (cut_case == 3) {
            if (tbox::Utilities::deq(areas[1],0.)) {
               flow_side[0] = 0;
            } else {
               flow_side[0] = 1;
            }
            if (tbox::Utilities::deq(areas[DIM],0.)) {
               flow_side[1] = 1;
            } else {
               flow_side[1] = 0;
            }
         }
         if (cut_case == 4) {
            if (tbox::Utilities::deq(areas[0],0.)) {
               flow_side[0] = 1;
            } else {
               flow_side[0] = 0;
            }
            if (tbox::Utilities::deq(areas[2*DIM-1],0.)) {
               flow_side[1] = 0;
            } else {
               flow_side[1] = 1;
            }
         }
         if (cut_case == 5) {
            if (tbox::Utilities::deq(areas[0],0.)) {
               flow_side[0] = 1;
            } else {
               flow_side[0] = 0;
            }
            if (tbox::Utilities::deq(areas[DIM],0.)) {
               flow_side[1] = 1;
            } else {
               flow_side[1] = 0;
            }
         }
         midpt[0] = (areas[DIM] + areas[2*DIM-1])/2.0;
         midpt[1] = (areas[0] + areas[1])/2.0;

      } // if DIM == 2

      if (DIM == 3) {
         TBOX_ERROR(d_object_name << ":calculateBoundaryNodeInformation()"
                          << "\nBoundary node computation does not work"
                            << "\nin 3D." << endl);
      }

      /*
       * Set centroid.
       */
      for (i = 0; i < DIM; i++) {
         if (cut_case == CUTCASE_UNDEFINED) {
            centroid[i] = tbox::IEEE::getSignalingNaN();
         } else {
            if (flow_side[i] == 0) {
               centroid[i] = lower[i] + midpt[i]*dx[i];
            } else if (flow_side[i] == 1) {
               centroid[i] = upper[i] - midpt[i]*dx[i];
            } else {
               TBOX_ERROR(d_object_name << ":calculateBoundaryNodeInformation()"
                          << "\nflow_side is not properly set." << endl);
            }
         }

         cut_cell.setFrontCentroid(centroid[i],i);

      }

      /*
       * Set shortest distance between boundary node and shape boundary
       */
      pdat::NodeIndex<DIM> ll(ic,hier::IntVector<DIM>(0,0));
      pdat::NodeIndex<DIM> lr(ic,hier::IntVector<DIM>(1,0));
      pdat::NodeIndex<DIM> ul(ic,hier::IntVector<DIM>(0,1));
      pdat::NodeIndex<DIM> ur(ic,hier::IntVector<DIM>(1,1));

      int num_bdry_nodes = cut_cell.getNumberBoundaryNodes();
      for (int j = 0; j < num_bdry_nodes; j++) {
         pdat::NodeIndex<NDIM> bdry_node = 
            cut_cell.getBoundaryNode(j).getIndex();

         double dist = EBGEOM_UNDEFINED;
         double distsq = EBGEOM_UNDEFINED;

         if (cut_case == 0) {
            if (bdry_node == ll || bdry_node == ul) {
               dist = dx[1] - areas[0];
            } else if (bdry_node == lr || bdry_node == ur) {
               dist = dx[1] - areas[1];
            } else {
               TBOX_ERROR(d_object_name << "::calculateBoundaryNodeInfo()"
                          << "\ndid not find boundary node in 'distance to "
                          << "shape' calculation" << endl);
            }
         } else if (cut_case == 1) {
            if (bdry_node == ll || bdry_node == lr) {
               dist = dx[0] - areas[DIM];
            } else if (bdry_node == ul || bdry_node == ur) {
               dist = dx[0] - areas[2*DIM-1];
            } else {
               TBOX_ERROR(d_object_name << "::calculateBoundaryNodeInfo()"
                          << "\ndid not find boundary node in 'distance to "
                          << "shape' calculation" << endl);
            }

         } else {

            /*
             * Set the distance to simply be the distance between 
             * the boundary node and the centroid.
             */
            if (bdry_node == ll) {
               node_loc[0] = lower[0];
               node_loc[1] = lower[1];
            } else if (bdry_node == lr) {
               node_loc[0] = upper[0];
               node_loc[1] = lower[1];
            } else if (bdry_node == ul) {
               node_loc[0] = lower[0];
               node_loc[1] = upper[1];
            } else if (bdry_node == ur) {
               node_loc[0] = upper[0];
               node_loc[1] = upper[1];
            } else {
               TBOX_ERROR(d_object_name << "::calculateBoundaryNodeInfo()"
                          << "\ndid not find boundary node in 'distance to "
                          << "shape' calculation" << endl);
            }
            
            distsq = 0.;
            for (i = 0; i < DIM; i++) {
               distsq += (node_loc[i] - centroid[i]) * 
                  (node_loc[i] - centroid[i]);
            } 
            dist = sqrt(distsq); 
         }

         BoundaryNode<DIM> bn1 =  cut_cell.getBoundaryNode(j);

         //cut_cell.getBoundaryNode(j).setDistanceToShapeBoundary(dist);
         bn1.setDistanceToBoundary(dist);
         cut_cell.setBoundaryNode(bn1,j);
         
         BoundaryNode<DIM> bn2 =  cut_cell.getBoundaryNode(j);
         
         dist = bn2.getDistanceToBoundary();

      }         
      
   } // if d_compute_areas_and_normal
}



/*
*************************************************************************
*                                                                       *
* Calculate the volume of the supplied cell.  This algorithm divides    *
* the cell into subcells using the maximum number of subdivides in each *
* direction.  Thus, each dimension is divided into 2^d_max_subdivides   *
* cells.  From these subcells, set an integer "inout" array centered    *
* at the nodes of the cell. Use the array of node flag values on the   *
* subcells to compute the volume fraction of the cell.                  *
*                                                                       *
* This method should return exactly the same result as the              *
* recursiveCalculateVolume() method because it uses the same algorithm  *
* implemented in a non-recursive form.  The user should be warned that  *
* the number of subcells can grow very quickly; for example, with a 3D  *
* problem with d_max_subdivides set to 3, each cell will be subdivided  *
* into 8x8x8 = 512 subcells.                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
EmbeddedBoundaryGeometry<DIM>::calculateVolume(
   const double *cell_lower,
   const double *cell_upper,
   double& error_estimate) const
{

   int i,k;
   double volume = 0.;

   /*
    * Compute the number of sub-cells (in each direction) from the max number 
    * of subdivides.  The total number of sub cells will be subcells^DIM.
    */
   int subcells = 1;
   for (i = 0; i < d_max_subdivides; i++) {
      subcells *= 2;
   }

   /*
    * Set the number of nodes, dx, and cell vol of each subcell
    */
   int nx[DIM];   
   double dx[DIM];
   int total_subnodes = 1;
   double subcell_vol = 1.;
   for (i = 0; i < DIM; i++) {
      nx[i] = subcells + 1;
      dx[i] = (cell_upper[i] - cell_lower[i])/(double)subcells;
      total_subnodes *= (subcells + 1);
      subcell_vol *= dx[i];
   }

   /*
    * Allocate space for node "inout" array, and set array by calling 
    * the shape's "isInside()" method.
    */
   int* node_flag = new int[total_subnodes];

   doNativeShapeInsideOutside(nx,dx,cell_lower,node_flag);
   
   /*
    * Compute total volume from subcell volumes using the following criteria:
    *   1. Compute subcell flag (flow, solid, cut) using:
    *           cell_flag = sum(node_flag) 
    *      where node_flag is the "flag" value at each corner node.
    *   2. If cell_flag = 2^DIM (solid)      volume += 0.
    *         cell_flag = 0      (flow)       volume += subcell_vol
    *         cell_flag > 0 && < 2^DIM (cut) volume += subcell_vol/2. 
    */
   int two_to_the_ndim = 1;
   for (i = 0; i < DIM; i++) {
      two_to_the_ndim *= 2;
   }
      
   int ijk, ijk_i, ijk_j, ijk_ij;
   int cell_flag = -1;
   int khi = 1;
   if (DIM == 3) {
      khi = nx[DIM-1]-1;
   }
   
   for (k = 0; k < khi; k++) {
      for (int j = 0; j < nx[1]-1; j++) {
         for (int i = 0; i < nx[0]-1; i++) {
            if (DIM == 2) {
               ijk = k*nx[1]*nx[0] + j*nx[0] + i;
               ijk_i = k*nx[1]*nx[0] + j*nx[0] + (i+1);
               ijk_j = k*nx[1]*nx[0] + (j+1)*nx[0] + i;
               ijk_ij = k*nx[1]*nx[0] + (j+1)*nx[0] + (i+1);
            
               cell_flag = 
                  node_flag[ijk] + node_flag[ijk_i] +
                  node_flag[ijk_j] + node_flag[ijk_ij];
            }
               
            if (DIM == 3) {
               ijk = (k+1)*nx[1]*nx[0] + j*nx[0] + i;
               ijk_i = (k+1)*nx[1]*nx[0] + j*nx[0] + (i+1);
               ijk_j = (k+1)*nx[1]*nx[0] + (j+1)*nx[0] + i;
               ijk_ij = (k+1)*nx[1]*nx[0] + (j+1)*nx[0] + (i+1);

               cell_flag += 
                  node_flag[ijk] + node_flag[ijk_i] +
                  node_flag[ijk_j] + node_flag[ijk_ij];
            }
            
            if (cell_flag == two_to_the_ndim) { // solid
               volume += 0.;
            } else if (cell_flag == 0) {  // flow
               volume += subcell_vol;
            } else {  // cut
               volume += 0.5 * subcell_vol;
               error_estimate += 0.25 * subcell_vol;
            }
         } // nx[0]
      } // nx[1]
   } // khi

   delete [] node_flag;
   
   return (volume);
}



/*
*************************************************************************
*                                                                       *
* Calculate the volume of the cut cell using a recursive algorithm.     *
* This involves recursively subdividing the cell and summing the        *
* subcell volumes until the minimum subcell volume (i.e. maximum        *
* number of subdivides) is reached.  The method returns the computed    *
* volume along with an estimate of the error.                           *
*                                                                       *
*************************************************************************
*/
template<int DIM> int  
EmbeddedBoundaryGeometry<DIM>::recursiveCalculateVolume(
   const double *cell_lower,
   const double *cell_upper,
   double& volume,
   double& error_estimate,
   int& subdivide_level) const
{
   int i,j;
   int return_val = 0;
   volume = 0.;

   /*
    * If the cell volume is less than the specified minimum subcell
    * volume, or if the maximum number of subdivides has been reached,
    * exit this method with volume = 0.
    */
   double cell_vol = 1.0;
   for (i = 0; i < DIM; i++) {
      cell_vol *= (cell_upper[i] - cell_lower[i]);
   }

   if (subdivide_level < d_max_subdivides) {
      
      /*
       * Form lower and upper extents of sub-boxes resulting from cell 
       * subdivision.
       *
       *   Bit shifting:  
       *    - Think of (1 << DIM) as 2**DIM 
       *      1d  1 << 1 = 100 = 2
       *      2d  1 << 2 = 010 = 4
       *      3d  1 << 3 = 001 = 8
       */      
      double lowers[(1 << DIM)][DIM];  // [ncorners][ndim]
      double uppers[(1 << DIM)][DIM];  // [ncorners][ndim]
      int classifications[(1 << DIM)];  // [ncorners]
      double midwy[DIM];
      
      subdivide_level++; // increment subdivide counter

      const int sub_boxes = (1 << DIM);
      
      int curr_subdivide_level = subdivide_level;

      for (i = 0; i < sub_boxes; i++) {
         for (j = 0; j < DIM; j++) {
            midwy[j] = 0.5*(cell_lower[j]+cell_upper[j]);
            lowers[i][j] = (i & (1 << j) ? cell_lower[j] : midwy[j]);
            uppers[i][j] = (i & (1 << j) ? midwy[j] : cell_upper[j]);
         }
         
         /*
          * Classify each of the sub boxes as FLOW, CUT, or SOLID.
          */
         classifications[i] = classifyCell(lowers[i], uppers[i]); 
         
      }
      
      /*
       * Add up the volumes based on classifications.  If it is CUT, 
       * call this method recursively.
       */
      double subvolume = cell_vol/(double)sub_boxes;
      
      for (i = 0; i < sub_boxes; i++) {
         if (classifications[i] == FLOW) {
            volume += subvolume;
         } else if (classifications[i] == CUT) {
            double cut_vol;
            int ret = recursiveCalculateVolume(lowers[i],
                                               uppers[i],
                                               cut_vol,
                                               error_estimate,
                                               subdivide_level);

            /*
             * If the returned sub-cell volume is zero, it means the sub
             * cell has reached our prescribed limit on size. Assume its 
             * split in half and update the error estimate.
             */
            double error_est = 0.;
            if (ret == 1) {
               cut_vol = 0.5 * subvolume;
               error_est = 0.25 * subvolume;
            } 
            volume += cut_vol;
            error_estimate += error_est;
         } else if (classifications[i] == SOLID) {
            volume += 0.;
         }

         /*
          * Reset subdivide_level back to original before moving to the
          * next sub-volume.  
          */

         subdivide_level = curr_subdivide_level;

      } // loop over sub-boxes
      
   } else {

      return_val = 1;

   }// subdivide_level < d_max_subdivides
 
   return (return_val);
}

/*
*************************************************************************
*                                                                       *
* Calculate the area lying within the computational domain of a cell.   *
* If ndim is 3d, then this function calculates the area of the faces of *
* the cell.  If ndim is 2d, it computes the length of the sides.  The   *
* faces are numbered from zero through ndim*2 and are of the form       *
* dim*2+side, where dim is the dimension and side is 0 (lower) or       *
* 1 (upper).                                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
EmbeddedBoundaryGeometry<DIM>::calculateArea(
   const double *cell_lower,
   const double *cell_upper,
   const int face,
   double& error_estimate) const
{
   int i;
   double area = 0.;

   double face_lower[DIM];
   double face_upper[DIM];   

   /*
    * Classify the face:
    *  face_dim - X, Y, or Z
    *  face_side - 0(lower) or 1(upper) face
    */
   int face_dim = EBGEOM_UNDEFINED;
   int face_side = EBGEOM_UNDEFINED;
   if (face == 0) {
      face_dim = 0;
      face_side = 0;
   } else if (face == 1) {
      face_dim = 0;
      face_side = 1;
   } else if (face == 2) {
      face_dim = 1;
      face_side = 0;
   } else if (face == 3) {
      face_dim = 1;
      face_side = 1;
   } else if (face == 4) {
      face_dim = 2;
      face_side = 0;
   } else if (face == 5) {
      face_dim = 2;
      face_side = 1;
   }
      
   double face_area = 1.0;
   for (i = 0; i < DIM; i++) {
      if (i == face_dim) {
         if (face_side == 0) {
            face_lower[i] = cell_lower[i];
            face_upper[i] = cell_lower[i];
         } else if (face_side == 1) {
            face_lower[i] = cell_upper[i];
            face_upper[i] = cell_upper[i];
         }
      } else {
         face_lower[i] = cell_lower[i];
         face_upper[i] = cell_upper[i];
         face_area *= face_upper[i] - face_lower[i];
      }
   }

   /*
    * Compute the number of sub-cells from the max number 
    * of subdivides.  The total number of sub cells will be subcells^DIM.
    */
   int subcells = 1;
   for (i = 0; i < d_max_subdivides; i++) {
      subcells *= 2;
   }

   /*
    * Set the number of nodes, dx, and cell vol of each subcell
    */
   int nx[DIM];   
   double dx[DIM];
   int total_subnodes = 1;
   double subface_area = 1.;
   for (i = 0; i < DIM; i++) {
      if (i == face_dim) {
         nx[i] = 1;
         dx[i] = 0.;
      } else {
         nx[i] = subcells + 1;
         dx[i] = (face_upper[i] - face_lower[i])/(double)subcells;
         total_subnodes *= (subcells + 1);
         subface_area *= dx[i];
      }
   }

   /*
    * Allocate space for node "inout" array, and set array by calling 
    * the shape's "isInside()" method.
    */
   int* node_flag = new int[total_subnodes];
   doNativeShapeInsideOutside(nx,dx,face_lower,node_flag);
   
   /*
    * Compute total area from subareas using the following criteria:
    *   1. Compute subarea flag (flow, solid, cut) using:
    *           cell_flag = sum(node_inout) 
    *      where node_flag is the "inout" value at each corner.
    *   2. If cell_flag = 2^DIM (solid)      area += 0.
    *         cell_flag = 0      (flow)      area += subcell_area
    *         cell_flag > 0 && < 2^DIM (cut) area += subcell_area/2. 
    */
   int two_to_the_ndim = 1;
   for (i = 0; i < DIM-1; i++) {
      two_to_the_ndim *= 2;
   }
      
   int ij, ij_i, ij_j, ij_ij;
   int cell_flag = -1;


   if (DIM == 2) {
      int face_nx[1];
      if (face_dim == 0) {
         face_nx[0] = nx[1];
      } else {
         face_nx[0] = nx[0];
      }
      for (int i = 0; i < face_nx[0]-1; i++) {
         cell_flag = node_flag[i] + node_flag[i+1];
            
         if (cell_flag == two_to_the_ndim) { // solid
            area += 0.;
         } else if (cell_flag == 0) {        // flow
            area += subface_area;
         } else {                            // cut
            area += 0.5 * subface_area;
            error_estimate += 0.25 * subface_area;
         }
         
      }
   } // DIM == 2

   if (DIM == 3) {
      int face_nx[2];
      if (face_dim == 0) {
         face_nx[0] = nx[1];
         face_nx[1] = nx[DIM-1];
      } else if (face_dim == 1) {
         face_nx[0] = nx[0];
         face_nx[1] = nx[DIM-1];
      } else if (face_dim == 2) {
         face_nx[0] = nx[0];
         face_nx[1] = nx[1];
      }
      for (int j = 0; j < face_nx[1]-1; j++) {
         for (int i = 0; i < face_nx[0]-1; i++) {
            ij    = j    *face_nx[0] + i;
            ij_i  = j    *face_nx[0] + (i+1);
            ij_j  = (j+1)*face_nx[0] + i;
            ij_ij = (j+1)*face_nx[0] + (i+1);
            
            cell_flag = 
               node_flag[ij] + node_flag[ij_i] +
               node_flag[ij_j] + node_flag[ij_ij];

            if (cell_flag == two_to_the_ndim) { // solid
               area += 0.;
            } else if (cell_flag == 0) {        // flow
               area += subface_area;
            } else {                            // cut
               area += 0.5 * subface_area;
               error_estimate += 0.25 * subface_area;
            }
         }
      }
   } // DIM == 3

   delete [] node_flag;

   return (area);


}


/*
*************************************************************************
*                                                                       *
* Calculate the area lying within the computational domain of a cell.   *
* If ndim is 3d, then this function calculates the area of the faces of *
* the cell.  If ndim is 2d, it computes the length of the sides.  The   *
* faces are numbered from zero through ndim*2 and are of the form       *
* dim*2+side, where dim is the dimension and side is 0 (lower) or       *
* 1 (upper).                                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::recursiveCalculateArea(
   const double *cell_lower,
   const double *cell_upper,
   const int face,
   double& area,
   double& error_estimate,
   int& subdivide_level) const
{
   int i,j;
   int return_val = 0;
   area = 0.;

   double face_lower[DIM];
   double face_upper[DIM];   

   /*
    * Classify the face:
    *  face_dim - X, Y, or Z
    *  face_side - 0(lower) or 1(upper) face
    */
   int face_dim = EBGEOM_UNDEFINED;
   int face_side = EBGEOM_UNDEFINED;
   if (face == 0) {
      face_dim = 0;
      face_side = 0;
   } else if (face == 1) {
      face_dim = 0;
      face_side = 1;
   } else if (face == 2) {
      face_dim = 1;
      face_side = 0;
   } else if (face == 3) {
      face_dim = 1;
      face_side = 1;
   } else if (face == 4) {
      face_dim = 2;
      face_side = 0;
   } else if (face == 5) {
      face_dim = 2;
      face_side = 1;
   }
      
   double face_area = 1.0;
   for (i = 0; i < DIM; i++) {
      if (i == face_dim) {
         if (face_side == 0) {
            face_lower[i] = cell_lower[i];
            face_upper[i] = cell_lower[i];
         } else if (face_side == 1) {
            face_lower[i] = cell_upper[i];
            face_upper[i] = cell_upper[i];
         }
      } else {
         face_lower[i] = cell_lower[i];
         face_upper[i] = cell_upper[i];
         face_area *= face_upper[i] - face_lower[i];
      }
   }

   if (subdivide_level < d_max_subdivides) {
      
      /*
       * Form lower and upper extents of sub-boxes resulting from face 
       * subdivision.
       *
       *   Bit shifting:  
       *    - Think of (1 << DIM) as 2**DIM 
       *      1d  1 << 1 = 100 = 2
       *      2d  1 << 2 = 010 = 4
       */      
      int area_ndim = DIM - 1;
      double lowers[(1 << area_ndim)][DIM];  // [ncorners][XYZ loc]
      double uppers[(1 << area_ndim)][DIM];  // [ncorners][XYZ loc]
      int classifications[(1 << area_ndim)];  // [ncorners]
      double midwy[DIM];
      
      subdivide_level++; // increment subdivide counter

      const int sub_boxes = (1 << area_ndim);
      
      int curr_subdivide_level = subdivide_level;

      for (i = 0; i < sub_boxes; i++) {
         for (j = 0; j < DIM; j++) {
            midwy[j] = 0.5*(face_lower[j] + face_upper[j]);
            lowers[i][j] = (i & (1 << j) ? face_lower[j] : midwy[j]);
            uppers[i][j] = (i & (1 << j) ? midwy[j] : face_upper[j]);
         }
         
         /*
          * Classify each of the sub boxes as FLOW, CUT, or SOLID.
          */
         classifications[i] = classifyCell(lowers[i], uppers[i]);
      }
      
      /*
       * Add up the areas based on classifications.  If it is CUT, 
       * call this method recursively.
       */
      double subarea = face_area/(double)sub_boxes;
      
      for (i = 0; i < sub_boxes; i++) {
         if (classifications[i] == FLOW) {
            area += subarea;
         } else if (classifications[i] == CUT) {
            double cut_area;
            int ret = recursiveCalculateArea(lowers[i],
                                             uppers[i],
                                             face,
                                             cut_area,
                                             error_estimate,
                                             subdivide_level);

            /*
             * If the returned sub-cell area is zero, it means the sub
             * cell has reached our prescribed limit on size. Assume its 
             * split in half and update the error estimate.
             */
            double error_est = 0.;
            if (ret == 1) {
               cut_area = 0.5 * subarea;
               error_est = 0.25 * subarea;
            } 
            area += cut_area;
            error_estimate += error_est;

         } else if (classifications[i] == SOLID) {
            area += 0.;
         }

         /*
          * Reset subdivide_level back to original before moving to the
          * next sub-volume.  
          */
         subdivide_level = curr_subdivide_level;

      } // loop over sub-boxes
      
   } else {

      return_val = 1;

   }

   return (return_val);
}

 
/*
*************************************************************************
*                                                                       *
* Given a cell lower and upper bounds, classify the cell as FLOW,       *
* SOLID, or CUT.  Note that this method may be applied to any           *
* dimension.  It is assumed that the "isInside()" function associated   *
* with each shape takes an ndim-vector and returns a boolean value.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> int 
EmbeddedBoundaryGeometry<DIM>::classifyCell(
   const double *lower, 
   const double *upper) const
{

   /*
    * Classify each of the corners of the cell and return the appropriate
    * flag depending on whether all corners are inside or outside of the
    * computational domain
    */

   bool classify_flow = true;
   bool classify_cut = false;
   bool classify_solid = false;   

   int corners = 1;

   for (int i=0; i<DIM; i++) {
      corners *=2;
   }
   const int ncorners = corners;
      
   bool is_inside = false;
   int i,c,n,inside_corner_ctr;

#ifdef USE_SINGLE_POINT_FOR_INOUT
   double corner[DIM];
#endif

   for (n = 0; n < d_shapes.getSize(); n++) {

#ifdef USE_ARRAY_FOR_INOUT
      int inout[ncorners];
      int nx[DIM];
      double dx[DIM];

      for (i = 0; i < DIM; i++) {
         nx[i] = 2;
         dx[i] = upper[i] - lower[i];
      }
      
      for (c = 0; c < ncorners; c++) {
         inout[c] = OUTSIDE;
      }
      
      d_shapes[n]->isInside(nx,
                            dx,
                            lower,
                            inout);

      inside_corner_ctr = 0;
      for (c = 0; c < ncorners; c++) {
         if (inout[c] == INSIDE) {
            is_inside = true;
         } else {
            is_inside = false;
         }
         if (is_inside) inside_corner_ctr++;
      }
#endif

#ifdef USE_SINGLE_POINT_FOR_INOUT
      inside_corner_ctr = 0;
      for (c = 0; c < ncorners; c++) {
         for (i = 0; i < DIM; i++) {
            corner[i] = (c & (1 << i) ? lower[i] : upper[i]);
         }
         is_inside = d_shapes[n]->isInside(corner);
         if (is_inside) inside_corner_ctr++;
      }
#endif
      if (inside_corner_ctr > ncorners) {
         TBOX_ERROR(d_object_name << ":classifyCell()"
                    << "\nSevere Error: the inside corner counter is "
                    << "\ngreater than the number of corners!" << endl);
      }

      /*
       * A cell with no corners inside will be flow; a cell with some
       * corners inside will be cut; and a cell with all corners inside
       * will be solid.
       *
       * Note that a cell may be flow for one shape, but cut or solid 
       * on another.  Or, it may be cut on one shape but solid on another.
       * For this reason, we establish the precedence:
       *
       * If cell is SOLID, this takes precedence over CUT and FLOW - 
       *    break out of the loop over shapes in this case.
       * If cell is CUT, this takes precedence over FLOW but not solid -
       *    we set classify_flow to false but do not reset solid and do
       *    not break out of the loop.
       * If cell is FLOW, it does not take precedence over either - we 
       *    do nothing in this case (classify_flow is set to true before
       *    the loop over shapes).
       */    
      if ((inside_corner_ctr > 0) && (inside_corner_ctr < ncorners)) {
         classify_flow = false;
         classify_cut = true;
      } else if (inside_corner_ctr == ncorners) {
         classify_flow = false;
         classify_cut = false;
         classify_solid = true;
         break;
      } 
      
   } // loop over shapes.

   int classification = -1;
   if (classify_cut) {
      classification = CUT;
   } else if (classify_flow) {
      classification = FLOW;      
   } else if (classify_solid) {
      classification = SOLID;
   } else {
      TBOX_ERROR(d_object_name << ":classifyCell()"
                 <<"\nCould not classify cell!!" << endl);
   }
   return(classification);

}


/*
*************************************************************************
*                                                                       *
*  Set the embedded boundary on the supplied level using data from      *
*  coarser level to initialize cell flags.  Note that this method       *
*  simply refines the cell flags to the finer level.  The steps to      *
*  step through the marked cells and compute cells that intersect       *
*  the triangulated surface are done via the "postprocessRefine()"      *
*  method, which gets called for each patch during refinement.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::refineEmbeddedBoundary(
   const tbox::Pointer<hier::PatchLevel<DIM> > level,
   const tbox::Pointer<hier::PatchLevel<DIM> > old_level,
   const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
#endif

   /*
    * If the "old_level" is null (i.e. not regridding) construct a schedule
    * that uses the coarser level data only.  If the "old_level" is non-null,
    * but no coarser level exists (e.g. hierarchy is null) construct a 
    * schedule that pulls data from the old level only. If the old_level is 
    * non-null and a coarser level exists, construct a schedule that pulls
    * data from both the old level and coarser level. 
    */ 
   int level_number = level->getLevelNumber();

   if (level_number == 0) {
      TBOX_ERROR(d_object_name << ":refineEmbeddedBoundary()\n"
                 << "Operating on level 0 - no coarser level available." 
                 << endl);
   }
     
   tbox::Pointer< xfer::RefineSchedule<DIM> > fill_sched;

   if (hierarchy.isNull()) {
      fill_sched = d_ebdry_refine_alg->createSchedule(level,
                                                      old_level,
                                                      this);
   } else {
      fill_sched = d_ebdry_refine_alg->createSchedule(level,
                                                      old_level,
                                                      level_number-1,
                                                      hierarchy,
                                                      this);
   }

   /*
    * Note: This step invokes the postprocessRefine() method
    */
   double time = 0.;
   fill_sched->fillData(time);

}

/*
*************************************************************************
*                                                                       *
* Set the embedded boundary on the physical boundaries (the boundary    *
* boxes) of the prescribed level.                                       *
*                                                                       *
* The implementation here takes advantage of the information stored     *
* in the BoundaryConditionsUtilities class.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::setEmbeddedBoundaryAtPhysicalBoundaries(
   const tbox::Pointer<hier::PatchLevel<DIM> > level)
{
   /*
    * Loop over patches on level and loop over boundary boxes.
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
         patch->getPatchGeometry();
      
      const hier::Box<DIM>& interior = patch->getBox();
      const hier::Index<DIM>& ifirst = interior.lower();
      const hier::Index<DIM>& ilast  = interior.upper();

      hier::IntVector<DIM> gcw_to_fill = d_ebdry_nghosts;

      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag = 
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);

      /*
       * These are the definitions for the different boundary conditions
       * defined in appu::CartesianBoundaryDefines.h.
       *
       * #define FLOW_BC            (90) 
       * #define REFLECT_BC         (91)
       * #define DIRICHLET_BC       (92)
       * #define NEUMANN_BC         (93)
       * 
       * The fortran will interpret the different "types" based on these 
       * values.
       */

      /*
       * Loop over NODE boundary boxes.
       */
      int btype = EBGEOM_UNDEFINED;
      if (DIM == 2) {
         btype = NODE2D_BDRY_TYPE;
      }            
      else if (DIM == 3) {
         btype = NODE3D_BDRY_TYPE;
      }            

      const tbox::Array<hier::BoundaryBox<DIM> >& node_bdry =
         pgeom->getCodimensionBoundary(btype);
      
      for (int i = 0; i < node_bdry.getSize(); i++) {
         
#ifdef DEBUG_CHECK_ASSERTIONS
         if (DIM == 2) {
            assert(node_bdry[i].getBoundaryType() == NODE2D_BDRY_TYPE);
         } else if (DIM == 3) {
            assert(node_bdry[i].getBoundaryType() == NODE3D_BDRY_TYPE);
         }
#endif
         

         int bnode_loc = node_bdry[i].getLocationIndex();
         
         hier::Box<DIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[i],
                                                  interior,
                                                  gcw_to_fill);
         const hier::Index<DIM> ibeg = fill_box.lower();
         const hier::Index<DIM> iend = fill_box.upper();         


         if (DIM == 2) {
            setebnode2d_(ifirst(0), ilast(0),
                         ifirst(1), ilast(1),
                         ibeg(0), iend(0),
                         ibeg(1), iend(1),
                         d_ebdry_nghosts(0), 
                         d_ebdry_nghosts(1),
                         bnode_loc,
                         d_node_bdry_cond[bnode_loc],
                         cell_flag->getPointer(),
                         cell_vol->getPointer());
         }
         

#if 0
         if (DIM == 3) {            
            setebnodebdry3d_(ifirst(0), ilast(0),
                             ifirst(1), ilast(1),
                             ifirst(2), ilast(2),
                             ibeg(0), iend(0),
                             ibeg(1), iend(1),
                             ibeg(2), iend(2),
                             d_ebdry_nghosts(0), 
                             d_ebdry_nghosts(1),
                             d_ebdry_nghosts(1),
                             bnode_loc,
                             d_node_bdry_conds[bnode_loc],
                             cell_flag->getPointer(),
                             cell_vol->getPointer());
         }
#endif

      }

      /*
       * Loop over EDGE boundary boxes.
       */

      if (DIM == 2) {
         btype = EDGE2D_BDRY_TYPE;
      }            
      else if (DIM == 3) {
         btype = EDGE3D_BDRY_TYPE;
      }            

      const tbox::Array<hier::BoundaryBox<DIM> >& edge_bdry =
         pgeom->getCodimensionBoundary(btype);

      for (int i = 0; i < edge_bdry.getSize(); i++) {
         
         hier::Box<DIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[i],
                                                             interior,
                                                  gcw_to_fill);
         const hier::Index<DIM> ibeg = fill_box.lower();
         const hier::Index<DIM> iend = fill_box.upper();

#if 0

         int bedge_loc = edge_bdry[i].getLocationIndex();
         

         if (DIM == 2) {
            
            setebedgebdry2d_(ifirst(0), ilast(0),
                             ifirst(1), ilast(1),
                             ibeg(0), iend(0),
                             ibeg(1), iend(1),
                             d_ebdry_nghosts(0), 
                             d_ebdry_nghosts(1),
                             bedge_loc,
                             d_edge_bdry_conds[bedge_loc],
                             cell_flag->getPointer(),
                             cell_vol->getPointer());
         }
         

         if (DIM == 3) {
            
            setebedgebdry3d_(ifirst(0), ilast(0),
                             ifirst(1), ilast(1),
                             ifirst(2), ilast(2),
                             ibeg(0), iend(0),
                             ibeg(1), iend(1),
                             ibeg(2), iend(2),
                             d_ebdry_nghosts(0), 
                             d_ebdry_nghosts(1),
                             d_ebdry_nghosts(1),
                             bedge_loc,
                             d_edge_bdry_conds[bedge_loc],
                             cell_flag->getPointer(),
                             cell_vol->getPointer());
         }
         
#endif

      }
         
      if (DIM == 3) {
         /*
          * Loop over FACE boundary boxes.
          */
         const tbox::Array<hier::BoundaryBox<DIM> >& face_bdry =
            pgeom->getCodimensionBoundary(FACE3D_BDRY_TYPE);
         
         for (int i = 0; i < face_bdry.getSize(); i++) {
            
         
         
#if 0

            int bface_loc = face_bdry[i].getLocationIndex();
            
            
            hier::Box<DIM> fill_box = pgeom->getBoundaryFillBox(face_bdry[i],
                                                                interior,
                                                                gcw_to_fill);
            
            const hier::Index<DIM> ibeg = fill_box.lower();
            const hier::Index<DIM> iend = fill_box.upper();
            
            setebfacebdry3d_(ifirst(0), ilast(0),
                             ifirst(1), ilast(1),
                             ifirst(2), ilast(2),
                             ibeg(0), iend(0),
                             ibeg(1), iend(1),
                             ibeg(2), iend(2),
                             d_ebdry_nghosts(0), 
                             d_ebdry_nghosts(1), 
                             d_ebdry_nghosts(2),
                             bface_loc,
                             d_face_bdry_conds[bface_loc],
                             cell_flag->getPointer(),
                             cell_vol->getPointer());
#endif
            
         } // loop over faces
      }

   } // loop over patches on level
}


/*
*************************************************************************
*                                                                       *
* Set parameters maintaining level ratio information for the hierarchy. *     
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::setLevelRatioInformation(
   const tbox::Pointer<hier::PatchLevel<DIM> > level)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(level.isNull()));
#endif

   /*
    *  Allocate level ratios and ratio to coarser arrays.
    */
   int level_num = level->getLevelNumber();
   
   if (level->inHierarchy() && level_num >= d_max_levels) {
      d_max_levels = level_num + 1;
      d_ratio_to_coarser.resizeArray(d_max_levels);
      d_level_ratios.resizeArray(d_max_levels);
   }

   /*
    * Set ratios for level 0.
    */
   d_ratio_to_coarser[0] = hier::IntVector<DIM>(1);
   d_level_ratios[0] = hier::IntVector<DIM>(1);

   /*
    * Set ratios for finer levels.
    */
   if (level_num > 0) {
      hier::IntVector<DIM> ratio = level->getRatioToCoarserLevel();
      d_ratio_to_coarser[level_num] = ratio;
      d_level_ratios[level_num] = 
         d_ratio_to_coarser[level_num] * d_level_ratios[level_num-1];
   }
   
}

/*
*************************************************************************
*                                                                       *
* Set node inside/outside array from native SAMRAI analytic shapes      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::doNativeShapeInsideOutside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* node_flag) const
{
   /*
    * Create temporary array to store node flag for each shape.
    */
   int total_nodes = 1;
   for (int i = 0; i < NDIM; i++) {
      total_nodes *= nx[i];
   }

   int* temp_node_flag = new int[total_nodes];

   /*
    * Loop over SAMRAI shapes and for each {
    *   - initialize all entries of temp_node_flag array to OUTSIDE
    *   - set temp_node_flag using shape's 'isInside()' method
    *   - copy temp_node_flag to supplied node_flag if temp_node_flag
    *     is INSIDE
    *   - clean up temporary array
    */  
   for (int nshape = 0; nshape < d_shapes.getSize(); nshape++) {      
      
      // initialize
      int n;
      for (n = 0; n < total_nodes; n++) {
         temp_node_flag[n] = OUTSIDE;
      }

      // set 
      d_shapes[nshape]->isInside(nx,
                                 dx,
                                 origin,
                                 temp_node_flag);
      
      // copy
      for (n = 0; n < total_nodes; n++) {
         if (temp_node_flag[n] == INSIDE) {
            node_flag[n] = INSIDE;
         }
      }

   } // loop over shapes

   // clean up
   delete [] temp_node_flag;
}


/*
*************************************************************************
*                                                                       *
* Write embedded boundary information - cell flag and cut cell          *
* information - to supplied file.                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::writeLevelEmbeddedBoundaryDataToFile(
   const tbox::Pointer<hier::PatchLevel<DIM> > level,
   const string& dirname) const
{

#ifdef HAVE_HDF5 

   /*
    * Generate directory to store mesh files.
    */
   bool write_to_dir = false;
   if (dirname.size() > 0) {
      write_to_dir = true;
      tbox::Utilities::recursiveMkdir(dirname);
   }
   
   /*
    * Form the filename.  It will be "ebmesh" appended
    * by the level number and processor id.
    *
    * i.e.   filename = "ebmesh-l<n>-p<pid>.hdf"
    */
   int i;
   int ln = level->getLevelNumber();
   int pid = tbox::MPI::getRank();
   tbox::pout << "\n  writing eb mesh to file = " << dirname 
        << "/ebmesh-l" << ln << "-p" << pid << "\n" << endl;
   
   string fileprefix = "ebmesh";
   const int size = fileprefix.length() + 16;
   char *buffer1 = new char[size];
   sprintf(buffer1, "%s-l%d-p%d", fileprefix.c_str(), ln, pid);
   string filename(buffer1);
   filename = filename + ".hdf";
   if (write_to_dir) filename = dirname + "/" + filename;
   
   /*
    * Open the HDF5 database with the supplied filename.
    */
   tbox::Pointer<tbox::HDFDatabase> db = new tbox::HDFDatabase("root");
   int stat = db->mount(filename,"W");
   if (stat < 0) {
     TBOX_ERROR(d_object_name << "writeLevelEmbeddedBoundaryDataToFile():" 
                << "\n Error opening HDF database: " << filename << endl);
   }

   /*
    * Write header information containing the box array for the 
    * level to verify matching when we read in the data.
    */          
   hier::BoxArray<DIM> patch_boxes = level->getBoxes();
   db->putDatabaseBoxArray("patch_boxes", patch_boxes);
   
   if (!d_grid_geometry.isNull()) {
      const double* domain_xlo  = d_grid_geometry->getXLower();
      const double* domain_xhi  = d_grid_geometry->getXUpper();
      db->putDoubleArray("domain_xlo",domain_xlo,DIM);
      db->putDoubleArray("domain_xhi",domain_xhi,DIM);
   }

   /*
    * Write embedded boundary data on the patches.  To keep the file size 
    * minimal, we only write the SOLID and CUT cells.  All
    * other cells are assumed to be FLOW.
    */
   
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());
      tbox::Pointer< pdat::IndexData<DIM,appu::CutCell<DIM> > > eboundary =
         patch->getPatchData(d_ebdry_data_id);
      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag =
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);
      
      tbox::List<hier::Index<DIM> > solid_cells;
      hier::Index<DIM> solid_cell(0);
      for (pdat::CellIterator<DIM> ic(patch->getBox()); ic; ic++) {
         if ((*cell_flag)(ic()) == SOLID) {
            for (i = 0; i < DIM; i++) {
               solid_cell(i) = ic()(i);  // cast CellIndex<DIM> to Index<DIM>
            }
            solid_cells.addItem(solid_cell);
         }
      }      
      
      /*
       * Create patch database to store solid and cut cell data.
       */
      char *buffer2 = new char[16];
      int patch_id = patch->getPatchNumber();
      sprintf(buffer2, "patch_db[%d]", patch_id);
      string name(buffer2);
      delete buffer2;
      tbox::Pointer<tbox::Database> patch_db = db->putDatabase(name);
      
      /*
       * Write solid cells for the patch.
       */
      int num_solid_cells = solid_cells.getNumberItems();
      patch_db->putInteger("num_solid_cells",num_solid_cells);
      
      if (num_solid_cells > 0) {
         
         // pack indices into an integer array
         tbox::Array<int> packed_indices;
         packed_indices.resizeArray(num_solid_cells*DIM);
         int cell = 0;
         for (typename tbox::List<hier::Index<DIM> >::Iterator 
                 li(solid_cells); li; li++) {
            for (i = 0; i < DIM; i++) {
               packed_indices[cell+i] = li()(i);
            }
            cell += DIM;
         }
         
         // write packed indices array to database
         patch_db->putIntegerArray("solid_cell_indices",packed_indices);
      }
      
      /*
       * Write cut cells for the patch
       */
      int num_cut_cells = eboundary->getNumberItems();
      patch_db->putInteger("num_cut_cells",num_cut_cells);
         
      int cut_cell_ctr = 0;

      for (typename pdat::IndexData<DIM,appu::CutCell<DIM> >::Iterator 
              bc(*eboundary); bc; bc++) {
         
         char *buffer3 = new char[16];
         sprintf(buffer3, "cut_cell[%d]", cut_cell_ctr);
         string name(buffer3);
         delete buffer3;
         tbox::Pointer<tbox::Database> cut_cell_db = 
            patch_db->putDatabase(name);
         
         tbox::Array<int> index;
         index.resizeArray(DIM);
         for (i = 0; i < DIM; i++) {
            index[i] = bc().getIndex()(i);
         }
         cut_cell_db->putIntegerArray("index",index);
         bc().putToDatabase(cut_cell_db);
         
         cut_cell_ctr++;

      }      

   } // loop over patches

   /*
    * Close the file.
    */
   db->unmount();

#endif
}

/*
*************************************************************************
*                                                                       *
* Read embedded boundary information - cell flag and cut cell           *
* information - from supplied file.                                     *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::readLevelEmbeddedBoundaryDataFromFile(
   const tbox::Pointer<hier::PatchLevel<DIM> > level,
   const string& dirname) const
{

#ifdef HAVE_HDF5
   /*
    * Form the filename.  It will be "ebmesh" appended
    * by the level number and processor id.
    *
    * i.e.   filename = "ebmesh-l<n>-p<pid>.hdf"
    */
   int ln = level->getLevelNumber();
   int pid = tbox::MPI::getRank();
   tbox::pout << "\n  reading eb mesh from file = " << dirname 
        << "/ebmesh-l" << ln << "-p" << pid << "\n" << endl;
   
   string fileprefix = "ebmesh";
   const int size = fileprefix.length() + 16;
   char *buffer1 = new char[size];
   sprintf(buffer1, "%s-l%d-p%d", fileprefix.c_str(), ln, pid);
   string filename(buffer1);
   filename = dirname + "/" + filename + ".hdf";
   
   /*
    * Open the HDF5 database with the supplied filename.
    */
   tbox::Pointer<tbox::HDFDatabase> db = new tbox::HDFDatabase(filename);
   if (db->mount(filename,"R") < 0) {
     TBOX_ERROR(d_object_name << "::readLevelEmbeddedBoundaryDataFromFile():"
                << "\n Error opening HDF database: " << filename << endl);
   }


   /*
    * Read patch boxes from file header.  This information is used to 
    * verify the written information and the read information is compatible.
    */
   hier::BoxArray<DIM> patch_boxes = db->getDatabaseBoxArray("patch_boxes");
   bool same = true;
   for (int p = 0; p < level->getNumberOfPatches(); p++) {
      if (p >= patch_boxes.getNumberOfBoxes()) {
         same = false;
         break;
      }
      hier::Box<DIM> rbox = patch_boxes(p);
      hier::Box<DIM> pbox = level->getPatch(p)->getBox();
      if (rbox != pbox) same = false;
   }
   if (!same) {
     TBOX_ERROR(d_object_name << "::readLevelEmbeddedBoundaryDataFromFile()"
                <<"\nBoxes in file: " << filename
                <<"\nare different than supplied Level " << ln
                <<" boxes. ***Exiting." << endl);
   }

   if (!d_grid_geometry.isNull()) {
      same = true;
      const double* domain_xlo  = d_grid_geometry->getXLower();
      const double* domain_xhi  = d_grid_geometry->getXUpper();
 
      double read_xlo[DIM];
      double read_xhi[DIM];
      db->getDoubleArray("domain_xlo",read_xlo,DIM);
      db->getDoubleArray("domain_xhi",read_xhi,DIM);

      for (int i = 0; i < DIM; i++) {
         if (!tbox::Utilities::deq(domain_xlo[i],read_xlo[i])) {
           TBOX_ERROR(d_object_name 
                      << "::readLevelEmbeddedBoundaryDataFromFile()"
                      <<"\nlevel xlo definition in: " << filename
                      <<"\nis different than supplied Level " << ln << " xlo"
                      <<"\n   level xlo[" << i << "]: " << domain_xlo[i]
                      <<"\n   read xlo[" << i << "]:  " << read_xlo[i]
                      << endl);
         }
         if (!tbox::Utilities::deq(domain_xhi[i],read_xhi[i])) {
           TBOX_ERROR(d_object_name 
                      << "::readLevelEmbeddedBoundaryDataFromFile()"
                      <<"\nlevel xhi definition in: " << filename
                      <<"\nis different than supplied Level " << ln << " xhi"
                      <<"\n   level xhi[" << i << "]: " << domain_xhi[i]
                      <<"\n   read xhi[" << i << "]:  " << read_xhi[i]
                      << endl);
         }
      }
   }
      
   /*
    * Read the embedded boundary information for each patch.  
    */   
   int i,d;
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      tbox::Pointer< pdat::IndexData<DIM,appu::CutCell<DIM> > > eboundary =
         patch->getPatchData(d_ebdry_data_id);
      tbox::Pointer< pdat::CellData<DIM,int> > cell_flag =
         patch->getPatchData(d_cell_flag_data_id);
      tbox::Pointer< pdat::CellData<DIM,double> > cell_vol =
         patch->getPatchData(d_cell_vol_data_id);

      /*
       * Determine name of patch database that holds solid and cut cell data.
       */
      char *buffer2 = new char[16];
      int patch_id = patch->getPatchNumber();
      sprintf(buffer2, "patch_db[%d]", patch_id);
      string name(buffer2);
      delete buffer2;
      tbox::Pointer<tbox::Database> patch_db = db->getDatabase(name);

      /*
       * Read solid cells from patch_db
       */
      int num_solid_cells = patch_db->getInteger("num_solid_cells");

      if (num_solid_cells > 0) {
         
         tbox::Array<int> packed_indices = 
            patch_db->getIntegerArray("solid_cell_indices");
         
         /*
          * Unpack solid cell indices and cell flag to SOLID for 
          * each index.
          */
         int cell_ctr = 0;
         for (i = 0; i < num_solid_cells; i++) {
            hier::Index<DIM> packed_index(0);
            for (d = 0; d < DIM; d++) {
               packed_index(d) = packed_indices[cell_ctr+d];
            }
            pdat::CellIndex<DIM> cell = packed_index;
            (*cell_flag)(cell) = SOLID;
            (*cell_vol)(cell) = 0.;
            cell_ctr += DIM;
         }
      }
      
         
      /*
       * Read cut cell info from patch_db
       */
      int num_cut_cells = patch_db->getInteger("num_cut_cells");
      for (i = 0; i < num_cut_cells; i++) {

         /*
          * Form cut cell database that will hold info about each cut cell
          */
         char *buffer3 = new char[16];
         sprintf(buffer3, "cut_cell[%d]", i);
         string name(buffer3);
         delete buffer3;
   
         /*
          * Access database and data in it.
          */
         tbox::Pointer<tbox::Database> cut_cell_db = 
            patch_db->getDatabase(name);

         tbox::Array<int> index_array = cut_cell_db->getIntegerArray("index");
         hier::Index<DIM> index(0);
         for (d = 0; d < DIM; d++) {
            index(d) = index_array[d];
         }

         /*
          * Create cut cell and add to eboundary
          */
         appu::CutCell<DIM> cc(index);
         cc.getFromDatabase(cut_cell_db);
         eboundary->appendItem(index, cc);

         /*
          * Set cell flag and volume fraction
          */
         pdat::CellIndex<DIM> cell = index;
         (*cell_flag)(cell) = CUT;
         (*cell_vol)(cell) = cc.getVolume();
         
      }
      
   } // loop over patches

   /*
    * Close the file.
    */
   db->unmount();

#endif
}


/*
*************************************************************************
*                                                                       *
* Allow user to set the grid geometry, if they did not supply it        *
* when the object was constructed.                                      *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::setGridGeometry(
   const tbox::Pointer<geom::CartesianGridGeometry<DIM> > grid_geom)
{
   d_grid_geometry = grid_geom;
}

/*
*************************************************************************
*                                                                       *
* Get data from input.                                                  *
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> ebdb,
   const bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(ebdb.isNull()));
#endif

   if (!is_from_restart) {
     d_flow_type = "EXTERNAL";
     if (ebdb->keyExists("flow_type")) {
       d_flow_type = ebdb->getString("flow_type");
     }
   }

   d_verbose = true;
   d_verbose = ebdb->getBoolWithDefault("verbose", d_verbose);
      
   d_read_ebdry = false;
   d_read_ebdry = ebdb->getBoolWithDefault("read_from_file", d_read_ebdry);
   
   d_write_ebdry = false;
   d_write_ebdry = ebdb->getBoolWithDefault("write_to_file", d_write_ebdry);
   
   if (ebdb->keyExists("dirname")) {
      d_ebdry_dirname = ebdb->getString("dirname");
   }
   
   int j;

   /*
    * Read in info about boundary definition.
    */


   if (ebdb->keyExists("CubesPatchInterface")) {

      /*
       * Boundary defined by Cubes.
       */
#ifndef HAVE_CUBES
      TBOX_ERROR("Cannot use Cubes - SAMRAI was not configured with it!");
#endif
      d_use_cubes = true;


   } else if (ebdb->keyExists("ElevenPatchInterface")) {

      /*
       * Boundary defined by Eleven.
       */
#ifndef HAVE_ELEVEN
      TBOX_ERROR("Cannot use Eleven - SAMRAI was not configured with it!");
#endif
      tbox::Pointer<tbox::Database> eleven_db = 
         ebdb->getDatabase("ElevenPatchInterface");
      d_use_eleven_inside_outside = 
         eleven_db->getBoolWithDefault("use_inside_outside",false);
      d_use_eleven_boundary_node = 
         eleven_db->getBoolWithDefault("use_boundary_node",false);
      if (!d_use_eleven_inside_outside && !d_use_eleven_boundary_node) {
         TBOX_ERROR(d_object_name << ": You must specify how you want to "
                    << "\nuse ELEVEN.  It will either compute nodal "
                    << "\ninside/outside info or compute Boundary nodes."
                    << "\nSet 'use_inside_outside' or 'use_boundary_node' in "
                    << "\nthe ElevenPatchInterface input to indicate which"
                    << "\nyou prefer." << endl);
      }
      if (d_use_eleven_inside_outside && d_use_eleven_boundary_node) {
         TBOX_ERROR(d_object_name << ": You have both 'use_inside_outside' and"
                    << "\n'use_boundary_node' specified in input.  You cannot "
                    << "\nuse both." << endl);
      }
      
   } else {
      
      /*
       * Boundary defined by native SAMRAI shapes.
       */
      if (ebdb->keyExists("Shapes")) {

         tbox::Pointer<tbox::Database> shapes_db = ebdb->getDatabase("Shapes");
   
         tbox::Array<string> shapes_keys = shapes_db->getAllKeys();
         int num_shapes = shapes_keys.getSize();
   
         if (d_flow_type == "INTERNAL" && num_shapes > 1) {
            TBOX_ERROR(d_object_name << ": Must supply only one shape when"
                       << " flow_type = INTERNAL" << endl);
         }
   
         /*
          * Each shape is expected to have input of the form:
          *   Shape {
          *     shape params
          *   }
          *
          * Here, we read in and construct each of the shapes provided
          * by the user.
          */
         d_shapes.resizeArray(num_shapes);
         for (j = 0; j < num_shapes; j++) {
            tbox::Pointer<tbox::Database> shape_db = 
               shapes_db->getDatabase(shapes_keys[j]);
            string type = shape_db->getString("type");
            tbox::Pointer< appu::EmbeddedBoundaryShape<DIM> > new_shape;
      
            if (type == "SPHERE") {
               string object_name = shapes_keys[j];
               new_shape =  new appu::EmbeddedBoundaryShapeSphere<DIM>(object_name, 
                                                                       shape_db);
            } else if (type == "POLYGON") {
               string object_name = shapes_keys[j];
               new_shape =  
                  new appu::EmbeddedBoundaryShapePolygon<DIM>(object_name, 
                                                              shape_db);
            } else {
               TBOX_ERROR("invalid shape type!!");
            }
            d_shapes[j] = new_shape;
      
         }
      
      }
   }
   
   /*
    * If SAMRAI is to be used for computing cutcell information, read
    * problem parameters controlling it here.
    */
   
   if (!d_use_cubes) {
      
      if (ebdb->keyExists("max_subdivides")) {
         d_max_subdivides = ebdb->getInteger("max_subdivides");
         if (d_verbose) {
            tbox::pout << "  d_max_subdivides = "
                       << d_max_subdivides << endl;
         }
      }
      
      d_compute_areas_and_normal = true;
      d_compute_areas_and_normal = 
         ebdb->getBoolWithDefault("compute_areas_and_normal",
                                  d_compute_areas_and_normal);
      
      d_use_recursive_algs = false;
      d_use_recursive_algs = 
         ebdb->getBoolWithDefault("use_recursive_algs",
                                  d_use_recursive_algs);
      
      d_compute_cutcell_index_data = true;
      d_compute_cutcell_index_data = 
         ebdb->getBoolWithDefault("compute_cutcell_index_data",
                                  d_compute_cutcell_index_data);
   }
   
      


   d_compute_boundary_node_data = false;
   d_compute_boundary_node_data = 
      ebdb->getBoolWithDefault("compute_boundary_node_data",
                               d_compute_boundary_node_data);

   /*
    * If boundary node data is to be used, specify that storage should
    * be allocated for it in the CutCell class.
    */
   if (d_compute_boundary_node_data) {
      CutCell<DIM>::enableBoundaryNodeStorage();
   }

   // Note: Currently the boundary node stuff only works in 2D.   
   if (DIM == 3 && d_compute_boundary_node_data) {
      TBOX_ERROR(d_object_name << ": Boundary nodes not working yet in 3D"
                 << "\nPlease set 'compute_boundary_node_data' in input"
                 << "to FALSE." << endl);
   }

   // Cannot compute boundary node data if compute_areas_and_normal is 
   // turned off.
   if (!d_compute_areas_and_normal && d_compute_boundary_node_data) {
      TBOX_ERROR(d_object_name << ": Cannot compute boundary node data"
                 << "\nif 'compute_areas_and_normal' is FALSE.  Please"
                 << "\nset 'compute_boundary_node_data' to FALSE." 
                 << endl);
   }



}




/*
*************************************************************************
*                                                                       *
* Read data from restart.                                               *
*                                                                       *
*************************************************************************
*/
   
template<int DIM> void 
EmbeddedBoundaryGeometry<DIM>::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> restart_db;
   if ( root_db->isDatabase(d_object_name) ) {
      restart_db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   int ver = restart_db->getInteger("APPU_EMBEDDED_BOUNDARY_GEOMETRY_VERSION");
   if (ver != APPU_EMBEDDED_BOUNDARY_GEOMETRY_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
                 << "Restart file version different than class version.");
   }

   d_flow_type = restart_db->getInteger("d_flow_type");
   tbox::Array<int> tmp_array = restart_db->getIntegerArray("d_ebdry_nghosts");
   for (int i = 0; i < DIM; i++) {
      d_ebdry_nghosts(i) = tmp_array[i];
   }

}

}
}
#endif
