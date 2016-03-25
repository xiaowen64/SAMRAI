//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/boundary/BoundaryDataTester.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Class to test usage of boundary utilities
//

#include "BoundaryDataTester.h"

#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"


//integer constants for boundary conditions
#define CHECK_BDRY_DATA  (1)
#include "CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA   (-9999)

// routines for managing boundary data
#if (NDIM == 2)
#include "CartesianBoundaryUtilities2.h"
typedef appu::CartesianBoundaryUtilities2 CartesianBoundaryUtilities;
#endif
#if (NDIM == 3)
#include "CartesianBoundaryUtilities3.h"
typedef appu::CartesianBoundaryUtilities3 CartesianBoundaryUtilities;
#endif

/*
*************************************************************************
*									*
* The constructor and destructor.                                       *
*									*
*************************************************************************
*/

BoundaryDataTester::BoundaryDataTester(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
#endif

   d_object_name = object_name;
   d_grid_geometry = grid_geom;

   d_variable_context = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("BOUNDARY_TEST");

   readVariableInputAndMakeVariables(input_db);

   setBoundaryDataDefaults();

   readBoundaryDataInput(input_db);

   postprocessBoundaryInput();

}

BoundaryDataTester::~BoundaryDataTester()
{
}

/*
*************************************************************************
*                                                                       *
* Set physical boundary values for each variable acording to input data.*
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double fill_time,
   const hier::IntVector<NDIM>& ghost_width_to_fill)
{
   tbox::plog << "\n\nFilling boundary data on patch = " << patch.getBox() << endl;
   tbox::plog << "ghost_width_to_fill = " << ghost_width_to_fill << endl;

   for (int iv = 0; iv < d_variables.getSize(); iv++) {

      tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = 
         patch.getPatchData(d_variables[iv], d_variable_context);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!cvdata.isNull());
#endif

      tbox::plog << "\n   iv = " << iv << " : " << d_variable_name[iv] << endl; 
      tbox::plog << "   depth = " << cvdata->getDepth() << endl;

      hier::IntVector<NDIM> fill_gcw = hier::IntVector<NDIM>::min(cvdata->getGhostCellWidth(),
                                          ghost_width_to_fill);

#if (NDIM == 3)
      CartesianBoundaryUtilities::
         fillFaceBoundaryData(d_variable_name[iv], cvdata,
                              patch,
                              fill_gcw,
                              ( (cvdata->getDepth() > 1)   ?
                                d_vector_bdry_face_conds   :
                                d_scalar_bdry_face_conds ),
                              d_variable_bc_values[iv]);
#endif  
      
      CartesianBoundaryUtilities::
         fillEdgeBoundaryData(d_variable_name[iv], cvdata,
                              patch,
                              fill_gcw,
                              ( (cvdata->getDepth() > 1)   ?
                                d_vector_bdry_edge_conds   :
                                d_scalar_bdry_edge_conds ),
                              d_variable_bc_values[iv]);

      CartesianBoundaryUtilities::
         fillNodeBoundaryData(d_variable_name[iv], cvdata,
                              patch,
                              fill_gcw,
                              ( (cvdata->getDepth() > 1)   ?
                                d_vector_bdry_node_conds   :
                                d_scalar_bdry_node_conds ),
                              d_variable_bc_values[iv]); 

   }

#if (NDIM == 2)
   checkBoundaryData(EDGE2D_BDRY_TYPE, patch, ghost_width_to_fill);
   checkBoundaryData(NODE2D_BDRY_TYPE, patch, ghost_width_to_fill);
#endif
#if (NDIM == 3)
   checkBoundaryData(FACE3D_BDRY_TYPE, patch, ghost_width_to_fill);
   checkBoundaryData(EDGE3D_BDRY_TYPE, patch, ghost_width_to_fill);
   checkBoundaryData(NODE3D_BDRY_TYPE, patch, ghost_width_to_fill);
#endif 

}

/*
*************************************************************************
*                                                                       *
* Set data for each variable on patch interior acording to input data.  *
*                                                                       *
*************************************************************************
*/
 
void BoundaryDataTester::initializeDataOnPatchInteriors(
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT(level_number == 0);
#endif

   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   level->allocatePatchData(d_patch_data_components);

   for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());

      for (int iv = 0; iv < d_variables.getSize(); iv++) {
         tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = 
            patch->getPatchData(d_variables[iv], d_variable_context);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!cvdata.isNull());
#endif
         for (int id = 0; id < cvdata->getDepth(); id++) {
            cvdata->fill(d_variable_interior_values[iv][id], 
                         patch->getBox(),
                         id);
         }
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Run boundary test:                                                    *
*                                                                       *
*  1) register boundary filling operation for each variable             *
*     with refine algorithm.                                            *
*                                                                       *
*  2) create communication schedule and fill data.                      *
*                                                                       *
*  3) check all patch boundary values for correctness.                  *
*                                                                       *
*************************************************************************
*/

int BoundaryDataTester::runBoundaryTest(
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT(level_number == 0);
#endif
   int d_fail_count = 0;

   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

   xfer::RefineAlgorithm<NDIM> boundary_fill;

   for (int iv = 0; iv < d_variables.getSize(); iv++) {
      int datid = 
         variable_db->mapVariableAndContextToIndex(d_variables[iv],
                                                   d_variable_context);

      boundary_fill.registerRefine(datid, datid, datid, NULL); 
   }

   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   boundary_fill.createSchedule(level, this)->fillData(0.0);

   return(d_fail_count);
}

/*
*************************************************************************
*                                                                       *
* Read variable data from input, create variables,                      *
* and map variables into variable database.                             *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::readVariableInputAndMakeVariables(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   tbox::Array<string> var_keys = db->getAllKeys();
   int nkeys = var_keys.getSize();

   int var_cnt = 0;
   for (int i = 0; i < nkeys; i++) {
      tbox::Pointer<tbox::Database> var_db = db->getDatabase(var_keys[i]);
      if (var_db->keyExists("name")) {
         var_cnt++;
      } 
   }

   d_variable_name.resizeArray(var_cnt);
   d_variable_depth.resizeArray(var_cnt);
   d_variable_num_ghosts.resizeArray(var_cnt);
   d_variable_interior_values.resizeArray(var_cnt);

   for (int i = 0; i < nkeys; i++) {

      tbox::Pointer<tbox::Database> var_db = db->getDatabase(var_keys[i]);

      if (var_keys[i] != "Boundary_data" && var_db->keyExists("name")) {

         if (var_db->keyExists("name")) {
            d_variable_name[i] = var_db->getString("name");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "Variable input error: No 'name' string found for "
               << "key = " << var_keys[i] << endl);
         }

         if (var_db->keyExists("depth")) {
            d_variable_depth[i] = var_db->getInteger("depth");
         } else {
            d_variable_depth[i] = 1; 
         }

         if (var_db->keyExists("num_ghosts")) {
            int* tmpg = d_variable_num_ghosts[i]; 
            var_db->getIntegerArray("num_ghosts", tmpg, NDIM);
         } else {
            d_variable_num_ghosts[i] = hier::IntVector<NDIM>(1); 
         }

         if (var_db->keyExists("interior_values")) {
            d_variable_interior_values[i].resizeArray(d_variable_depth[i]); 
            var_db->getDoubleArray("interior_values", 
                                   d_variable_interior_values[i].getPointer(),
                                   d_variable_depth[i]);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "Variable input error: No 'interior_values' entry found for "
               << "key = " << var_keys[i] << endl);
         }

      }

   }

   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

   d_variables.resizeArray(d_variable_name.getSize());

   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {
      d_variables[iv] = new pdat::CellVariable<NDIM,double>(d_variable_name[iv],
                                                 d_variable_depth[iv]);

      int datid = 
         variable_db->registerVariableAndContext(d_variables[iv],
                                                 d_variable_context,
                                                 d_variable_num_ghosts[iv]);

      d_patch_data_components.setFlag(datid);
   }

}

/*
*************************************************************************
*                                                                       *
* Set all boundary data to bogus default values for error checking.     *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::setBoundaryDataDefaults()
{
   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

#if (NDIM == 2)
   d_master_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_node_bdry_edge.resizeArray(NUM_2D_NODES);

   for (int ni = 0; ni < NUM_2D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
   }
#endif

#if (NDIM == 3)
   d_master_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_scalar_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_vector_bdry_face_conds.resizeArray(NUM_3D_FACES);
   for (int fi = 0; fi < NUM_3D_FACES; fi++) {
      d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_edge_bdry_face.resizeArray(NUM_3D_EDGES);
   for (int ei = 0; ei < NUM_3D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_node_bdry_face.resizeArray(NUM_3D_NODES);

   for (int ni = 0; ni < NUM_3D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
   }
#endif

   d_variable_bc_values.resizeArray(d_variable_name.getSize());
   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {
#if (NDIM == 2)
      d_variable_bc_values[iv].resizeArray(NUM_2D_EDGES*d_variable_depth[iv]);
#endif
#if (NDIM == 3)
      d_variable_bc_values[iv].resizeArray(NUM_3D_FACES*d_variable_depth[iv]);
#endif
      tbox::MathUtilities<double>::setArrayToSignalingNaN(d_variable_bc_values[iv]);
   }

}

/*
*************************************************************************
*                                                                       *
* Functions to read boundary information from input database.           *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::readDirichletBoundaryDataEntry(
   tbox::Pointer<tbox::Database> db, 
   string& db_name, 
   int bdry_location_index)
{
   readBoundaryDataStateEntry(db, db_name, bdry_location_index);
}

void BoundaryDataTester::readNeumannBoundaryDataEntry(
   tbox::Pointer<tbox::Database> db,
   string& db_name,
   int bdry_location_index)
{
   readBoundaryDataStateEntry(db, db_name, bdry_location_index);
}

void BoundaryDataTester::readBoundaryDataStateEntry(
   tbox::Pointer<tbox::Database> db,
   string& db_name,
   int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(d_variable_bc_values.getSize() == d_variable_name.getSize());
#endif

   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {

#ifdef DEBUG_CHECK_ASSERTIONS
#if (NDIM == 2)
      TBOX_ASSERT(d_variable_bc_values[iv].getSize() == 
             NUM_2D_EDGES*d_variable_depth[iv]);
#endif
#if (NDIM == 3)
      TBOX_ASSERT(d_variable_bc_values[iv].getSize() == 
             NUM_3D_FACES*d_variable_depth[iv]);
#endif
#endif

       if (db->keyExists(d_variable_name[iv])) {
          int depth = d_variable_depth[iv];
          tbox::Array<double> tmp_val(0);
          tmp_val = db->getDoubleArray(d_variable_name[iv]);
          if (tmp_val.getSize() < depth) {
             TBOX_ERROR(d_object_name << ": "
                << "Insufficient number of " 
                << d_variable_name[iv] << " values given in " 
                << db_name << " input database." << endl);
          }
          for (int id = 0; id < depth; id++) {
             d_variable_bc_values[iv][bdry_location_index*depth+id] = 
                tmp_val[id];
          }
       } else {
          TBOX_ERROR(d_object_name << ": "
             << d_variable_name[iv] 
             << " entry missing from " << db_name
             << " input database. " << endl);
       }

   }

}
  
void BoundaryDataTester::readBoundaryDataInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   hier::IntVector<NDIM> periodic = d_grid_geometry->getPeriodicShift();
   int num_per_dirs = 0;
   for (int id = 0; id < NDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < NDIM) {

      if (db->keyExists("Boundary_data")) {

         tbox::Pointer<tbox::Database> bdry_db = db->getDatabase("Boundary_data");

         CartesianBoundaryUtilities::
            readBoundaryInput(this,
                              bdry_db,
#if (NDIM == 3)
                              d_master_bdry_face_conds,
#endif
                              d_master_bdry_edge_conds,
                              d_master_bdry_node_conds,
                              periodic);

      } else {
         TBOX_ERROR(d_object_name << ": "
            << "Key data 'Boundary_data' not found in input. " << endl);
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Postprocess boundary data from input values                           *
* to make setting and checking easier.                                  *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::postprocessBoundaryInput()
{
#if (NDIM == 2)
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == REFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = FLOW_BC;
      }
   }
   for (int i = 0; i < NUM_2D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_edge[i] =
            appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }
#endif
#if (NDIM == 3)
   for (int i = 0; i < NUM_3D_FACES; i++) {
      d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
      d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];

      if (d_master_bdry_face_conds[i] == REFLECT_BC) {
         d_scalar_bdry_face_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_3D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
         d_edge_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                                            i, d_master_bdry_edge_conds[i]);
      }
   }

   for (int i = 0; i < NUM_3D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }
#endif

}

/*
*************************************************************************
*                                                                       *
* Check boundary values on patch for correctness.                       *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::checkBoundaryData(
   int btype, 
   const hier::Patch<NDIM>& patch, 
   const hier::IntVector<NDIM>& ghost_width_to_check)
{
#ifdef DEBUG_CHECK_ASSERTIONS
#if (NDIM == 2) 
   TBOX_ASSERT(btype == EDGE2D_BDRY_TYPE ||
          btype == NODE2D_BDRY_TYPE);
#endif
#if (NDIM == 3)
   TBOX_ASSERT(btype == FACE3D_BDRY_TYPE ||
          btype == EDGE3D_BDRY_TYPE ||
          btype == NODE3D_BDRY_TYPE);
#endif
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const tbox::Array<hier::BoundaryBox<NDIM> > bdry_boxes =
      pgeom->getCodimensionBoundary(btype);

   for (int i = 0; i < bdry_boxes.getSize(); i++ ) {
      hier::BoundaryBox<NDIM> bbox = bdry_boxes[i];
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
#endif
      int bloc = bbox.getLocationIndex();

      for (int iv = 0; iv < d_variables.getSize(); iv++) {
         tbox::Pointer< pdat::CellData<NDIM,double> > cvdata =
            patch.getPatchData(d_variables[iv], d_variable_context);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!cvdata.isNull());
#endif

         int depth = d_variable_depth[iv];

         int bscalarcase, bvectorcase, refbdryloc;
#if (NDIM == 2) 
         if (btype == EDGE2D_BDRY_TYPE) {
            bscalarcase = d_scalar_bdry_edge_conds[bloc];
            bvectorcase = d_vector_bdry_edge_conds[bloc];
            refbdryloc = bloc;
         } else { // btype == NODE2D_BDRY_TYPE
            bscalarcase = d_scalar_bdry_node_conds[bloc];
            bvectorcase = d_vector_bdry_node_conds[bloc];
            refbdryloc = d_node_bdry_edge[bloc];
         }
#endif
#if (NDIM == 3)
         if (btype == FACE3D_BDRY_TYPE) {
            bscalarcase = d_scalar_bdry_face_conds[bloc];
            bvectorcase = d_vector_bdry_face_conds[bloc];
            refbdryloc = bloc;
         } else if (btype == EDGE3D_BDRY_TYPE) {
            bscalarcase = d_scalar_bdry_edge_conds[bloc];
            bvectorcase = d_vector_bdry_edge_conds[bloc];
            refbdryloc = d_edge_bdry_face[bloc];
         } else { // btype == NODE3D_BDRY_TYPE
            bscalarcase = d_scalar_bdry_node_conds[bloc];
            bvectorcase = d_vector_bdry_node_conds[bloc];
            refbdryloc = d_node_bdry_face[bloc];
         }
#endif

         int data_id = hier::VariableDatabase<NDIM>::getDatabase()->
            mapVariableAndContextToIndex(d_variables[iv], d_variable_context);

         int num_bad_values = 0;

         if (depth == 1) {

            num_bad_values =
            CartesianBoundaryUtilities::
               checkBdryData(d_variable_name[iv],
                             patch,
                             data_id,
                             0,
                             ghost_width_to_check,
                             bbox,
                             bscalarcase,
                             d_variable_bc_values[iv][refbdryloc]);
#if (TESTING == 1)
            if (num_bad_values > 0) {
               d_fail_count++; 
               tbox::perr << "\nBoundary Test FAILED: \n"
                    << "     " << num_bad_values << " bad "
                    << d_variable_name[iv] << " values found for"
                    << "     boundary type " << btype << " at location " 
                                                      << bloc << endl;
            }
#endif

         } else {
            for (int id = 0; id < depth; id++) {
               int vbcase = bscalarcase;
#if (NDIM == 2)
               if (btype == EDGE2D_BDRY_TYPE) {
                  if ( (id == 0 && (bloc == XLO || bloc == XHI)) ||
                       (id == 1 && (bloc == YLO || bloc == YHI)) ) {
                     vbcase = bvectorcase;
                  }
               } else {
                   if ( (id == 0 && bvectorcase == XREFLECT_BC) ||
                        (id == 1 && bvectorcase == YREFLECT_BC) ) {
                      vbcase = bvectorcase;
                   }
               }
#endif
#if (NDIM == 3)
               if (btype == FACE3D_BDRY_TYPE) {
                  if ( (id == 0 && (bloc == XLO || bloc == XHI)) ||
                       (id == 1 && (bloc == YLO || bloc == YHI)) ||
                       (id == 2 && (bloc == ZLO || bloc == ZHI)) ) {
                     vbcase = bvectorcase;
                  }
               } else {
                   if ( (id == 0 && bvectorcase == XREFLECT_BC) ||
                        (id == 1 && bvectorcase == YREFLECT_BC) ||
                        (id == 2 && bvectorcase == ZREFLECT_BC) ) {
                      vbcase = bvectorcase;
                   }
               }
#endif

               num_bad_values =
               CartesianBoundaryUtilities::
                  checkBdryData(d_variable_name[iv],
                                patch,
                                data_id,
                                id,
                                ghost_width_to_check,
                                bbox,
                                vbcase,
                                d_variable_bc_values[iv][refbdryloc*depth+id]);
#if (TESTING == 1)
               if (num_bad_values > 0) {
                  d_fail_count++; 
                  tbox::perr << "\nBoundary Test FAILED: \n"
                       << "     " << num_bad_values << " bad "
                       << d_variable_name[iv] << " values found for"
                       << "     boundary type " << btype << " at location " 
                                                         << bloc << endl;
               }
#endif
            
            }  // for (int id = 0; id < depth; id++)

         }  // else

      }   // for (int iv = 0; iv < d_variables.getSize(); iv++)

   }  // for (int i = 0; i < bdry_boxes.getSize(); i++ )

}

/*
*************************************************************************
*                                                                       *
* Write all class data members to specified output stream.              *
*                                                                       *
*************************************************************************
*/

void BoundaryDataTester::printClassData(ostream &os) const
{
   int i,j;
   os << "\nBoundaryDataTester::printClassData..." << endl;
   os << "BoundaryDataTester: this = " << (BoundaryDataTester*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::CartesianGridGeometry<NDIM>*)d_grid_geometry << endl;

   if (!d_variable_context.isNull()) {
      os << "d_variable_context = "
         << d_variable_context->getName() << endl;
   } else {
      os << "d_variable_context = NULL" << endl;
   }

   os << "\nVariables ...\n" << endl;
   for (i = 0; i < d_variable_name.getSize(); i++) {
      os << "Variable " << i << endl;  
      os << "   name       = " << d_variable_name[i] << endl;  
      os << "   depth      = " << d_variable_depth[i] << endl;  
      os << "   num_ghosts = " << d_variable_num_ghosts[i] << endl;  
      os << "   interior_values = " << d_variable_interior_values[i][0]; 
      for (j = 1; j < d_variable_depth[i]; j++) { 
         os << " ,  " << d_variable_interior_values[i][j];
      }
      os << endl;
   }

   os << "\n   Boundary condition data... " << endl;

#if (NDIM == 2)
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      if (d_master_bdry_edge_conds[j] == DIRICHLET_BC ||
          d_master_bdry_edge_conds[j] == NEUMANN_BC) {
         for (i = 0; i < d_variable_name.getSize(); i++) {
            os << d_variable_name[i] << " bdry edge value[" << j << "] = "
               << d_variable_bc_values[i][j*d_variable_depth[i]];
            for (int id = 1; id < d_variable_depth[i]; id++) {
               os << " , " 
                  << d_variable_bc_values[i][j*d_variable_depth[i]+id];
            }
            os << endl;
         }
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_edge[" << j << "] = "
         << d_node_bdry_edge[j] << endl;
   }
#endif
#if (NDIM == 3)
   for (j = 0; j < d_master_bdry_face_conds.getSize(); j++) {
      os << "\n       d_master_bdry_face_conds[" << j << "] = "
         << d_master_bdry_face_conds[j] << endl;
      os << "       d_scalar_bdry_face_conds[" << j << "] = "
         << d_scalar_bdry_face_conds[j] << endl;
      os << "       d_vector_bdry_face_conds[" << j << "] = "
         << d_vector_bdry_face_conds[j] << endl;
      if (d_master_bdry_face_conds[j] == DIRICHLET_BC) {
         for (i = 0; i < d_variable_name.getSize(); i++) {
            os << d_variable_name[i] << " bdry edge value[" << j << "] = "
               << d_variable_bc_values[i][j*d_variable_depth[i]];
            for (int id = 1; id < d_variable_depth[i]; id++) {
               os << " , "
                  << d_variable_bc_values[i][j*d_variable_depth[i]+id];
            }
            os << endl;
         }
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      os << "       d_edge_bdry_face[" << j << "] = "
         << d_edge_bdry_face[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_face[" << j << "] = "
         << d_node_bdry_face[j] << endl;
   }
#endif


}
