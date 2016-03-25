//
// File:        VisItDataWriter.C
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:    $Revision: 727 $
// Modified:    $Date: 2005-11-10 17:20:10 -0800 (Thu, 10 Nov 2005) $
// Description: Writes data files for visualization by VisIt
//


#ifndef included_appu_VisItDataWriter_C
#define included_appu_VisItDataWriter_C

#include "VisItDataWriter.h"

#ifdef HAVE_HDF5

#include "tbox/Array.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "MultiblockPatchHierarchy.h"
#include "NodeData.h"
#include "NodeDataFactory.h"
#include "tbox/MPI.h"
#include "tbox/IEEE.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#define VISIT_NAME_BUFSIZE (128)
#define VISIT_UNDEFINED_INDEX (-1)

// for parallel runs, VISIT_MASTER writes single summary file 
// with information from all processors.
#define VISIT_MASTER  (0)   


// used by VisIt to track version of VisIt Data Writer
#define VISIT_DATAWRITER_VERSION_NUMBER 2.0


extern "C" {
  void cpfdat2buf3d_(
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  float*,  double*, const int&);
  void cpddat2buf3d_(
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  double*,  double*, const int&);
  void cpidat2buf3d_(
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  const int& , const int& , const int& , 
  int*,  double*, const int&);
}
extern "C" {
  void cpfdat2buf2d_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  float*,  double*, const int&);
  void cpddat2buf2d_(
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& , 
  const int& , const int& , 
  double*,  double*, const int&);
  void cpidat2buf2d_(
  const int& , const int& , 
  const int& , const int& ,
  const int& , const int& , 
  const int& , const int& ,
  int*,  double*, const int&);
}

namespace SAMRAI {
    namespace appu {

template<int DIM> bool VisItDataWriter<DIM>::s_summary_file_opened = false;

   
/*
*************************************************************************
*                                                                       *
* The constructor --- sets default object state.                        *
*                                                                       *
*************************************************************************
*/

template<int DIM>  VisItDataWriter<DIM>::VisItDataWriter(
   const string& object_name,
   const string& dump_directory_name,
   int number_procs_per_file,
   bool is_multiblock)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(number_procs_per_file > 0);
#endif

   if ((DIM < 2) || (DIM > 3)) {
      TBOX_ERROR(
	 "VisItDataWriter<DIM>::VisItDataWriter<DIM>"
	 << "\n          VisItDataWriter only works for"
	 << "\n          2D or 3D data" << endl);
   }

   d_object_name = object_name;

   d_default_derived_writer = NULL;
   d_materials_writer = NULL;

   d_number_working_slaves = VISIT_UNDEFINED_INDEX;
   d_file_cluster_size = number_procs_per_file;
   d_number_file_clusters = VISIT_UNDEFINED_INDEX;
   d_my_file_cluster_number = VISIT_UNDEFINED_INDEX;
   d_file_cluster_leader = false;
   d_my_rank_in_file_cluster = VISIT_UNDEFINED_INDEX;
   d_number_files_this_file_cluster = VISIT_UNDEFINED_INDEX;

   d_scaling_ratios.resizeArray(1);
   d_scaling_ratios[0] = hier::IntVector<DIM>(1);

   d_number_visit_variables = 0; 
   d_number_visit_variables_plus_depth = 0; 
   d_number_species = 0;

   d_time_step_number = VISIT_UNDEFINED_INDEX;
   d_grid_type = VISIT_CARTESIAN;
   d_top_level_directory_name = dump_directory_name;
   d_summary_filename = "summary.samrai";
   d_number_levels = 1;

   d_worker_min_max = (patchMinMaxStruct*)NULL;

   d_is_multiblock = is_multiblock;

}

/*
*************************************************************************
*                                                                       *
* The destructor implicitly deallocates the list of plot data items.    *
*                                                                       *
*************************************************************************
*/

template<int DIM>  VisItDataWriter<DIM>::~VisItDataWriter()
{
   /*
    * De-allocate min/max structs for each variable.
    */
   if (d_worker_min_max != (patchMinMaxStruct*)NULL) 
      delete [] d_worker_min_max;

   
   for (typename tbox::List<VisItItem>::Iterator ipi(d_plot_items); 
        ipi; ipi++) {
      for (int comp = 0; comp < VISIT_MAX_NUMBER_COMPONENTS; comp++) {
         if (ipi().d_master_min_max[comp] != (patchMinMaxStruct*)NULL) 
            delete [] ipi().d_master_min_max[comp];
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set default derived data writer.                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::setDefaultDerivedDataWriter(
   VisDerivedDataStrategy<DIM>* derived_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(derived_writer != (VisDerivedDataStrategy<DIM>*)NULL);
#endif
   d_default_derived_writer = derived_writer;
}

/*
*************************************************************************
*                                                                       *
* Set materials data writer object.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::setMaterialsDataWriter(
   VisMaterialsDataStrategy<DIM>* materials_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(materials_writer != (VisMaterialsDataStrategy<DIM>*)NULL);
#endif
   d_materials_writer = materials_writer;
}

/*
*************************************************************************
*                                                                       *
* Register (non-derived) plot quantities: scalar, vector or tensor.     *
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::registerPlotQuantity(
   const string& variable_name, 
   const string& variable_type, 
   const int patch_data_index,
   const int start_depth_index,
   const double scale_factor,
   const string& variable_centering)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable_name.empty());
   assert(!variable_type.empty());
   assert(patch_data_index >= -1);
   assert(start_depth_index >= 0);
#endif

   /*
    * Check for name conflicts with existing registered variables.
    */
   for (typename tbox::List<VisItItem>::Iterator
          ipi(d_plot_items); ipi; ipi++) {
      if (ipi().d_var_name == variable_name) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity()"
            << "\n    Attempting to register variable with name " 
            << variable_name << "\n    more than once." << endl);
      }
   }

   /*
    * Create a plot item and initialize its characteristics.
    */
   VisItItem plotitem;

   initializePlotItem(plotitem, 
                      variable_name,
                      variable_type,
                      patch_data_index,
                      start_depth_index,
                      scale_factor,
                      variable_centering);

   d_number_visit_variables++;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.appendItem(plotitem);
   
}

/*
*************************************************************************
*                                                                       *
* Register derived plot quantities: scalar, vector, or tensor.  If no   *
* derived data strategy is specified and no default derived data        *
* strategy is set, an error will result.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::registerDerivedPlotQuantity(
   const string& variable_name,
   const string& variable_type,
   VisDerivedDataStrategy<DIM>* derived_writer,
   double scale_factor,
   const string& variable_centering)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable_name.empty());
   assert(!variable_type.empty());
#endif

   /*
    * Check for name conflicts with existing registered variables.
    */
   for (typename tbox::List<VisItItem>::Iterator
          ipi(d_plot_items); ipi; ipi++) {
      if (ipi().d_var_name == variable_name) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerDerivedPlotQuantity()"
            << "\n    Attempting to register variable with name " 
            << variable_name << "\n    more than once." << endl);
      }
   }

   /*
    * Create a plot item and initialize its characteristics.
    */
   VisItItem plotitem;

   /*
    * Derived data is packed by the user.  Thus, we specify a dummy
    * patch data id here.
    */
   int patch_data_index = VISIT_UNDEFINED_INDEX;  
   int start_depth_index = 0;
   initializePlotItem(plotitem, 
                      variable_name,
                      variable_type,
                      patch_data_index,
                      start_depth_index,
                      scale_factor,
                      variable_centering);

   /*
    * Set characteristics for derived variable.
    */
   plotitem.d_is_derived = true;

   if (derived_writer == NULL) {
      if (d_default_derived_writer == NULL) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerDerivedPlotQuantity"
                    << "\n    no derived data writer specified for variable:"
                    <<        variable_name
                    << "\n    and no default derived data writer set." 
                    << endl);
      } else {
         plotitem.d_derived_writer = d_default_derived_writer;
      }
   } else {
      plotitem.d_derived_writer = derived_writer;
   }      

   d_number_visit_variables++;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.appendItem(plotitem);
}

/*
*************************************************************************
*                                                                       *
* Reset previously-registered scalar/vector variable to new data id     *
* and depth index on the given level.  This allows the use of different * 
* patch data ids for the same quantity on different hierarchy levels.   *
* We check to make sure that the factory at the given index is          *
* defined and consistent with the original registration.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::resetLevelPlotQuantity(
   const string& variable_name,
   const int level_number,
   const int patch_data_index,
   const int start_depth_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable_name.empty());
   assert(level_number >= 0);
   assert(patch_data_index >= -1);
   assert(start_depth_index >= 0);
#endif

   /*
    * Verify the supplied patch data index has the same type and centering
    * as the plot item its replacing.
    */
   tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
      hier::VariableDatabase<DIM>::getDatabase()->
                           getPatchDescriptor()->
                              getPatchDataFactory(patch_data_index);

   bool found_type = false;
   variable_data_type vdt = VISIT_DATA_TYPE_BAD;
   variable_centering vc = VISIT_CENTERING_BAD;
#ifdef HAVE_FLOAT
   if (!found_type) {
      tbox::Pointer< pdat::CellDataFactory<DIM,float> > ffactory = factory;
      if (!ffactory.isNull()) {
         vdt = VISIT_FLOAT;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      tbox::Pointer< pdat::NodeDataFactory<DIM,float> > ffactory = factory;
      if (!ffactory.isNull()) {
         vdt = VISIT_FLOAT;
         vc = VISIT_NODE;
         found_type = true;
      }
   }
#endif
   if (!found_type) {
      tbox::Pointer< pdat::CellDataFactory<DIM,double> > dfactory = factory;
      if (!dfactory.isNull()) {
         vdt = VISIT_DOUBLE;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      tbox::Pointer< pdat::NodeDataFactory<DIM,double> > dfactory = factory;
      if (!dfactory.isNull()) {
         vdt = VISIT_DOUBLE;
         vc = VISIT_NODE;
         found_type = true;
      }
   }
   if (!found_type) {
      tbox::Pointer< pdat::CellDataFactory<DIM,int> > ifactory = factory;
      if (!ifactory.isNull()) {
         vdt = VISIT_INT;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      tbox::Pointer< pdat::NodeDataFactory<DIM,int> > ifactory = factory;
      if (!ifactory.isNull()) {
         vdt = VISIT_INT;
         vc = VISIT_NODE;
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter<DIM>::resetLevelPlotQuantity()"
                 << "\n     Unable to determine type and centering"
                 << "\n     of supplied patch data index."
                 << "\n     ***Exiting" << endl);
   }

   /*
    * Find variable in the list of maintained vars.
    */
   bool found_var = false;
   for (typename tbox::List<VisItItem>::Iterator 
        ipi(d_plot_items); (!found_var && ipi); ipi++) {
      if (ipi().d_var_name == variable_name) {

         /*
          * Check to make sure supplied variable has same type and
          * centering as registered one.
          */
         if ((ipi().d_var_data_type != vdt) || 
             (ipi().d_var_centering != vc)) {
            TBOX_ERROR("VisItDataWriter<DIM>::resetLevelPlotQuantity()"
                       << "\n     The supplied patch data id has a different"
                       << "\n     type and centering from the one originally"
                       << "\n     registered.  hier::Variable name: "
                       << variable_name
                       << "\n     ***Exiting" << endl);
         }
         if (level_number >= ipi().d_level_patch_data_index.getSize()) {
            ipi().d_level_patch_data_index.resizeArray(level_number+1);
         }
         ipi().d_level_patch_data_index[level_number] = patch_data_index;
         ipi().d_start_depth_index = start_depth_index;
      }
   }

   if (!found_var) {
      TBOX_ERROR("VisItDataWriter<DIM>::resetLevelPlotQuantity()"
                 << "\n     Could not find the variable: "
                 << variable_name
                 << "\n     in the list of registered plot items."
                 << "\n     Be sure registerPlotQuantity() has been"
                 << "\n     called for the variable."
                 << "\n     ***Exiting" << endl);
   }
      
}

/*
*************************************************************************
*                                                                       *
* Register node coordinates of deformed (moving) grids.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::registerNodeCoordinates(
   const int patch_data_index,
   const int start_depth_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(patch_data_index >= -1);
   assert(start_depth_index >= 0);
#endif

   /*
    * Check to make sure "Coords" variable has not already been registered.
    */
   for (typename tbox::List<VisItItem>::Iterator
          ipi(d_plot_items); ipi; ipi++) {
      if (ipi().d_var_name == "Coords") {
         TBOX_ERROR("VisItDataWriter<DIM>::registerNodeCoordinates()"
            << "\n   Coordinates registered more than once." << endl);
      }
   }

   /*
    * Set the grid type for the visit data
    */
   d_grid_type = VISIT_DEFORMED;
   
   
   /*
    * Verify the supplied patch data index is a valid NODE-centered
    * float or double and has a depth of at least DIM
    */
   tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
      hier::VariableDatabase<DIM>::getDatabase()->
                           getPatchDescriptor()->
                              getPatchDataFactory(patch_data_index);

   bool found_type = false;
   int var_depth = VISIT_UNDEFINED_INDEX;
   if (!found_type) {
#ifdef HAVE_FLOAT
      tbox::Pointer< pdat::NodeDataFactory<DIM,float> > ffactory = factory;
      if (!ffactory.isNull()) {
         var_depth = ffactory->getDefaultDepth();
         found_type = true;
      }
#endif
   }
   if (!found_type) {
      tbox::Pointer< pdat::NodeDataFactory<DIM,double> > dfactory = factory;
      if (!dfactory.isNull()) {
         var_depth = dfactory->getDefaultDepth();
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerNodeCoordinates"
                 << "\n     This variable is NOT a node centered"
                 << "\n     float or double type, which is required."
                 << "\n     ***Exiting" << endl);
   }

   int end_depth = start_depth_index + DIM;
   if (var_depth < (end_depth)) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerNodeCoordinates"
                 << "\n     This variable has depth: " << var_depth
                 << "\n     It must be a VECTOR type and therefore"
                 << "\n     have depth at least DIM + start_depth_index = " 
                 << end_depth
                 << "\n     ***Exiting" << endl);
   }
   

   /*
    * Create the coords plot item.
    */
   VisItItem plotitem;
   
   string var_name = "Coords";
   string var_type = "VECTOR";
   double scale_factor = 1.0;
   string var_cent = "NODE";
   
   initializePlotItem(plotitem, 
                      var_name,
                      var_type,
                      patch_data_index,
                      start_depth_index,
                      scale_factor,
                      var_cent);

   plotitem.d_is_deformed_coords = true;
   
   /*
    * We need to reset the variable name, because it has to be written with
    * a special form to the VisIt readible HDF file.
    */
   char temp_buf[VISIT_NAME_BUFSIZE];
   for (int i = 0; i < plotitem.d_depth; i++) {
      sprintf(temp_buf, ".%02d",i);
      plotitem.d_visit_var_name[i] = var_name + temp_buf;
   }
      
   
   /*
    * In this method, we assume the scale factor is always 1.0.  If a 
    * user would like to choose a scale factor, use the 
    * "registerSingleNodeCoordinate()" method.
    */
   plotitem.d_coord_scale_factor.resizeArray(DIM);
   for (int i = 0; i < DIM; i++) {
      plotitem.d_coord_scale_factor[i] = 1.0;
   }

   d_number_visit_variables++;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.appendItem(plotitem);

}

/*
*************************************************************************
*                                                                       *
* Register node coordinates of deformed (moving) grids.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::registerSingleNodeCoordinate(
   const int coordinate_number,
   const int patch_data_index,
   const int depth_index,
   const double scale_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(coordinate_number >= 0 && coordinate_number < DIM);
   assert(patch_data_index >= -1);
   assert(depth_index >= 0);
#endif

   /*
    * Set the grid type for the visit data
    */
   d_grid_type = VISIT_DEFORMED;
   
   
   /*
    * Verify the supplied patch data index is a valid NODE-centered
    * float or double
    */
   tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
      hier::VariableDatabase<DIM>::getDatabase()->
                           getPatchDescriptor()->
                              getPatchDataFactory(patch_data_index);

   bool found_type = false;
   if (!found_type) {
#ifdef HAVE_FLOAT
      tbox::Pointer< pdat::NodeDataFactory<DIM,float> > ffactory = factory;
      if (!ffactory.isNull()) {
         found_type = true;
      }
#endif
   }
   if (!found_type) {
      tbox::Pointer< pdat::NodeDataFactory<DIM,double> > dfactory = factory;
      if (!dfactory.isNull()) {
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerSingleNodeCoordinate"
                 << "\n     This variable is NOT a node centered"
                 << "\n     float or double type, which is required."
                 << "\n     ***Exiting" << endl);
   }

   /*
    * Create the coords plot item.  If its the first time this method
    * is called, initialize the "Coords" plot variable.  If it has 
    * already been called before (i.e. coordinate_number > 0) then just
    * reset plot variable parameters as necessary.
    */
   if (coordinate_number == 0) {

      /*
       * Check to make sure "Coords" variable has not already been registered.
       */
      for (typename tbox::List<VisItItem>::Iterator
              ipi(d_plot_items); ipi; ipi++) {
         if (ipi().d_var_name == "Coords") {
            TBOX_ERROR("VisItDataWriter<DIM>::registerSingleNodeCoordinate()"
                       << "\n   Coordinate registered more than once." 
                       << endl);
         }
      }

      VisItItem plotitem;
      
      string var_name = "Coords";
      string var_type = "SCALAR";
      string var_cent = "NODE";
   
      initializePlotItem(plotitem, 
                         var_name,
                         var_type,
                         patch_data_index,
                         depth_index,
                         scale_factor,
                         var_cent);

      plotitem.d_is_deformed_coords = true;

      /*
       * We need to reset the variable name, because it has to be written with
       * a special form to the VisIt readible HDF file.
       */
      char temp_buf[VISIT_NAME_BUFSIZE];
      for (int i = 0; i < plotitem.d_depth; i++) {
         sprintf(temp_buf, ".%02d",i);
         plotitem.d_visit_var_name[i] = var_name + temp_buf;
      }

      plotitem.d_coord_scale_factor.resizeArray(DIM);
      plotitem.d_coord_scale_factor[coordinate_number] = scale_factor;  
      d_number_visit_variables++;
      d_number_visit_variables_plus_depth += plotitem.d_depth;    

      d_plot_items.appendItem(plotitem);

   } else {
      
      for (typename tbox::List<VisItItem>::Iterator
              ipi(d_plot_items); ipi; ipi++) {

         if (ipi().d_is_deformed_coords) {

            ipi().d_var_type = VISIT_VECTOR;
            ipi().d_depth = DIM;
            ipi().d_visit_var_name.resizeArray(DIM);
            
            string var_name = "Coords";
            char temp_buf[VISIT_NAME_BUFSIZE];
            sprintf(temp_buf, ".%02d",coordinate_number);
            ipi().d_visit_var_name[coordinate_number] = var_name + temp_buf;

            ipi().d_coord_scale_factor[coordinate_number] = scale_factor;
            d_number_visit_variables_plus_depth += 1;
         }
      }
   }

}
   
/*
*************************************************************************
*                                                                       *
* Register material names -- names of all the materials (not species)   *
* being used in the application.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::registerMaterialNames(
   const tbox::Array<string>& material_names)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(material_names.getSize() > 0);
#endif
   /*
    * Check if we have already tried to register materials.
    */
   if (d_materials_names.getSize() > 0) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerMaterialNames"
                 << "\n    This method has been called more than once."
                 << "\n    The material names may not change during the"
                 << "\n    simulation.  ***Exiting" << endl);
   }


   /*
    * Register each of the material names as a plot item with material
    * characteristics.
    */
   int num_materials = material_names.getSize();
   d_materials_names.resizeArray(num_materials);
   for (int i = 0; i < num_materials; i++) {
      if (material_names[i].empty()) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerMaterialNames"
                    << "\n    Material: "  << i
                    << "\n    has an empty name.  A name must be supplied."
                    << "\n    ***Exiting" << endl);
      }      

      d_materials_names[i] = material_names[i];
      
      VisItItem plotitem;

      /*
       * We impose the criteria that materials must be scalar cell centered
       * double quantities.  The user must supply a method to pack
       * this information.  See header for more info.
       */
      string var_type = "SCALAR";
      int patch_data_index = VISIT_UNDEFINED_INDEX;
      int start_depth_index = 0;
      double scale_factor = 1.0;
      string var_cent = "CELL";

      initializePlotItem(plotitem, 
                         material_names[i],
                         var_type,
                         patch_data_index,
                         start_depth_index,
                         scale_factor,
                         var_cent);

      plotitem.d_isa_material = true;
      plotitem.d_material_name = material_names[i];
 
      d_plot_items.appendItem(plotitem);
   }
}

/*
*************************************************************************
*                                                                       *
* Register species names for a given material.                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::registerSpeciesNames(
   const string& material_name,
   const tbox::Array<string>& species_names) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!material_name.empty());
   assert(species_names.getSize() > 0);
#endif

   /*
    * Be sure we have already registered materials.
    */
   if (d_materials_names.getSize() == 0) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerSpeciesNames"
                 << "\n    No materials have yet been registered."
                 << "\n    Be sure the 'registerMaterialNames()'"
                 << "\n    is called before this method.  ***Exiting" 
                 << endl);
   }


   bool found_material = false;
   for (typename tbox::List<VisItItem>::Iterator ipi(d_plot_items); ipi; ipi++) {
      if (ipi().d_material_name == material_name) {
         found_material = true;
      }
   }
   if (!found_material) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerSpeciesNames"
                 << "\n    material name = " << material_name
                 << "\n    has not been registered." << endl);
   }

   /*
    * Find the material in the list of plot items
    */
   VisItItem *material_item = (VisItItem*)NULL;
   for (typename tbox::List<VisItItem>::Iterator ipi(d_plot_items); ipi; ipi++) {
      if ((ipi().d_material_name == material_name) &&
                             ipi().d_isa_material) {
         if (ipi().d_species_names.getSize() > 0) {
            TBOX_ERROR("VisItDataWriter<DIM>::registerSpeciesNames"
                       << "\n    material name = " << material_name
                       << "\n    registerSpeciesNames has already been"
                       << "\n    called for this material.  It is not"
                       << "\n    possible to reset species during the"
                       << "\n    simulation.  ***Exiting" << endl);
         } 
         material_item = &(ipi());
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(material_item != (VisItItem*)NULL);
#endif

   d_number_species += species_names.getSize();
   material_item->d_species_names = species_names;

   /*
    * Create a plot variable for each species of the material.
    */
   for (int i = 0; i < species_names.getSize(); i++) {

      if (species_names[i].empty()) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerSpeciesNames"
                    << "\n    species i = " << i
                    << "\n    for material: " << material_item->d_material_name
                    << "\n    is empty.  ***Exiting" << endl);
      }

      /*
       * Species are associated with a material.  Hence, we give them
       * the same characteristics (centering, type, etc) as material
       * variables.
       */
      VisItItem plotitem;

      string var_type = "SCALAR";
      int patch_data_index = VISIT_UNDEFINED_INDEX;
      int start_depth_index = 0;
      double scale_factor = 1.0;
      string var_cent = "CELL";

      initializePlotItem(plotitem,
                         species_names[i],
                         var_type,
                         patch_data_index,
                         start_depth_index,
                         scale_factor,
                         var_cent);      

      plotitem.d_isa_species = true;

      plotitem.d_species_name = species_names[i];
      plotitem.d_parent_material_pointer = material_item;
      
      d_plot_items.appendItem(plotitem);

   }

}



/*
*************************************************************************
*                                                                       *
* Private method which initializes a VisIt variable based on user       *
* input.
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::initializePlotItem(
   VisItItem& plotitem,
   const string& variable_name, 
   const string& variable_type, 
   const int patch_data_index,
   const int start_depth_index,
   const double scale_factor,
   const string& variable_centering)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable_name.empty());
   assert(!variable_type.empty());
   assert(patch_data_index >= -1);
   assert(start_depth_index >= 0);
#endif

   plotitem.d_var_name = variable_name;

   /*
    * Set variable type.
    */
   if (variable_type == "SCALAR") {
      plotitem.d_var_type = VISIT_SCALAR;
      plotitem.d_depth = 1;
   } else if (variable_type == "VECTOR") {
      plotitem.d_var_type = VISIT_VECTOR;
      plotitem.d_depth = DIM;
   } else if (variable_type == "TENSOR") {
      plotitem.d_var_type = VISIT_TENSOR;
      plotitem.d_depth = DIM*DIM;
   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity"
                 << "\n    variable_type " << variable_type
                 << "\n    is unsupported.  You must use SCALAR, VECTOR, or" 
                 << "\n    TENSOR.  Exiting***" << endl);
   }

   /*
    * Check to make sure we have not exceeded max allowed components.
    */
   int num_old_components = d_number_visit_variables_plus_depth
      + d_materials_names.getSize() // number materials 
      + d_number_species;
   
   int new_num_components = num_old_components + plotitem.d_depth;
   if (new_num_components > VISIT_MAX_NUMBER_COMPONENTS) {
      TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity"
                 << "\n     Unable to register this quantity because it"
                 << "\n     the maximum number of variables allowed in"
                 << "\n     the VisItWriter was reached:"
                 << "\n       current num variables:" 
                 << num_old_components
                 << "\n       variable depth: " 
                 << plotitem.d_depth
                 << "\n     MAX_NUMBER_COMPONENTS: "
                 << VISIT_MAX_NUMBER_COMPONENTS
                 << "\n     Contact SAMRAI team for assistance." << endl);
   }

   /*
    * Set variable centering.  If the variable supplied a valid patch
    * index, we check the factory of the variable first to try to 
    * determine centering from that.  If we cannot determine
    * the type, we use the "variable_centering" provided by the user.
    * Note that we only set the centering for variables that are not
    * derived, materials, or species, since these variables do not have
    * a supplied patch data index because they are packed by the user.
    */
   bool found_type = false;
   int var_depth = 0;
   if (patch_data_index >= 0) {
      
      tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
         hier::VariableDatabase<DIM>::getDatabase()->
         getPatchDescriptor()->
         getPatchDataFactory(patch_data_index);
      if (factory.isNull()) {
         TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity"
                    << "\n    patch data array index = " << patch_data_index 
                    << "\n    for variable = " << variable_name
                    << "\n    is invalid" << endl);
      } else {

         if (!found_type) {
#ifdef HAVE_FLOAT
            tbox::Pointer< pdat::CellDataFactory<DIM,float> > ffactory = factory;
            if (!ffactory.isNull()) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_FLOAT;
               var_depth = ffactory->getDefaultDepth();
               found_type = true;
            }
#endif
         }
         if (!found_type) {
            tbox::Pointer< pdat::CellDataFactory<DIM,double> > dfactory = factory;
            if (!dfactory.isNull()) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_DOUBLE;
               var_depth = dfactory->getDefaultDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            tbox::Pointer< pdat::CellDataFactory<DIM,int> > ifactory = factory;
            if (!ifactory.isNull()) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_INT;
               var_depth = ifactory->getDefaultDepth();
               found_type = true;
            }
         }
         if (!found_type) {
#ifdef HAVE_FLOAT
            tbox::Pointer< pdat::NodeDataFactory<DIM,float> > ffactory = factory;
            if (!ffactory.isNull()) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_FLOAT;
               var_depth = ffactory->getDefaultDepth();
               found_type = true;
            }
#endif
         }
         if (!found_type) {
            tbox::Pointer< pdat::NodeDataFactory<DIM,double> > dfactory = factory;
            if (!dfactory.isNull()) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_DOUBLE;
               var_depth = dfactory->getDefaultDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            tbox::Pointer< pdat::NodeDataFactory<DIM,int> > ifactory = factory;
            if (!ifactory.isNull()) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_INT;
               var_depth = ifactory->getDefaultDepth();
               found_type = true;
            }
         }

         /*
          * Make sure variable depth is sufficient for the specified type
          * (SCALAR/VECTOR/TENSOR) with the start depth index.
          */
         int end_depth = start_depth_index + plotitem.d_depth;
         if (var_depth < end_depth) {
            TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity"
                       << "\n    The variable: " << variable_name
                       << "\n    has insufficient depth for the type" 
                       << "\n    and start depth index registered."
                       << "\n      var_type:           " << variable_type
                       << "\n      start_depth_index:  " << start_depth_index
                       << "\n      required min depth: " << end_depth
                       << endl);
         }         

      } // factory.isNull
 
   } // valid patch data index
      
   if (!found_type) {
      if (variable_centering == "CELL") {
         plotitem.d_var_centering = VISIT_UNKNOWN_CELL;
      } else if (variable_centering == "NODE") {
         plotitem.d_var_centering = VISIT_UNKNOWN_NODE;
      } else {
         TBOX_ERROR("VisItDataWriter<DIM>::registerPlotQuantity"
                    << "\n     Unable to determine the centering for"
                    << "\n     this variable; it must be supplied."
                    << "\n     The variable_centering argument is " 
                    << variable_centering 
                    << "\n     Possible entries are CELL or NODE" 
                    << "\n     ***Exiting" << endl);
      }
      
      /*
       * When the var centering is specified by the user, we assume
       * it is of type double.  This implies the user should NOT try
       * to register their undefined variables with type float or int. 
       */
      plotitem.d_var_data_type = VISIT_DOUBLE;
   }

   /*
    * Set the patch data index.
    */
   plotitem.d_patch_data_index = patch_data_index;
   plotitem.d_level_patch_data_index.resizeArray(d_number_levels);
   for (int ln = 0; ln < d_number_levels; ln++) {
      plotitem.d_level_patch_data_index[ln] = patch_data_index;
   }

   plotitem.d_visit_var_name.resizeArray(plotitem.d_depth);
   char temp_buf[VISIT_NAME_BUFSIZE];
   int i;
   for (i = 0; i < plotitem.d_depth; i++) {
      if (plotitem.d_depth == 1) {
         plotitem.d_visit_var_name[i] = variable_name;
      } else {
         sprintf(temp_buf, ".%02d",i);
         plotitem.d_visit_var_name[i] = variable_name + temp_buf;
      }
   }
   

   plotitem.d_scale_factor = scale_factor;
   plotitem.d_start_depth_index = start_depth_index;

   /*
    * Initialize min/max information.
    */
   for (i = 0; i < VISIT_MAX_NUMBER_COMPONENTS; i++) {
      plotitem.d_master_min_max[i] = (patchMinMaxStruct*)NULL;
   }
   
   /*
    * Set derived, coords, material, and species information NULL.  
    * If the variable is any of these types, this information
    * should be set by the appropriate registration functions.
    */
   plotitem.d_is_derived = false;
   plotitem.d_derived_writer = (VisDerivedDataStrategy<DIM>*)NULL;
   plotitem.d_is_deformed_coords = false;

   plotitem.d_isa_material = false;
   plotitem.d_materials_writer = (VisMaterialsDataStrategy<DIM>*)NULL;
   
   plotitem.d_isa_species = false;
   plotitem.d_parent_material_pointer = (VisItItem*)NULL;

}




/*
*************************************************************************
*                                                                       *
* Private functions for parallel runs which serve as barriers to enable *
* orderly writing of cluster files by passing a baton from the current  *
* writer to the next proc in cluster, and so on.                        *
*                                                                       *
*************************************************************************
*/

#define VISIT_FILE_CLUSTER_WRITE_BATON 17

template<int DIM> void VisItDataWriter<DIM>::dumpWriteBarrierBegin()
{
   int x[1], proc_before_me, len = 1;

   if (d_file_cluster_leader) {
      return;
   } else {
      proc_before_me = (d_my_file_cluster_number * d_file_cluster_size) + 
                        d_my_rank_in_file_cluster - 1;

      tbox::MPI::recv((int *)x, 
                     len, 
                     proc_before_me, 
                     false, 
                     VISIT_FILE_CLUSTER_WRITE_BATON);
   }
}

template<int DIM> void VisItDataWriter<DIM>::dumpWriteBarrierEnd()
{
   int x[1], proc_after_me;
   int num_procs = tbox::MPI::getNodes();
   proc_after_me = (d_my_file_cluster_number * d_file_cluster_size) +
                    d_my_rank_in_file_cluster + 1;
   if (proc_after_me < num_procs) {
      tbox::MPI::send(x, 1, proc_after_me, false, 
                     VISIT_FILE_CLUSTER_WRITE_BATON);
   }
}

/*
*************************************************************************
*                                                                       *
* Write plot data from given hierarchy to HDF file                      * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::writePlotData(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   int time_step_number,
   double simulation_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(time_step_number >= 0);
   assert(!d_top_level_directory_name.empty());
#endif

   if (time_step_number <= d_time_step_number) {
      TBOX_ERROR("VisItDataWriter<DIM>::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n    time step number: " << time_step_number
         << " is <= last time step number: " << d_time_step_number
         << endl);
   }
   d_time_step_number = time_step_number;

   if ((d_materials_names.getSize() > 0) && 
       (d_materials_writer == NULL)) {
      TBOX_ERROR("VisItDataWriter<DIM>::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n    setMaterialsDataWriter() has not been called," 
         << "\n    this method must be called when using materials."
         << endl);
   }

   if (d_is_multiblock) {
      d_number_levels = 0;
      int finest_level_num = hierarchy->getFinestLevelNumber();
      for (int ln = 0; ln <= finest_level_num; ln++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
            hierarchy->getPatchLevel(ln);
         for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               d_number_levels++;
            }
         }
      }
   } else {

      d_number_levels = hierarchy->getNumberLevels();

   }
   if (d_number_levels > d_scaling_ratios.getSize()) {
      d_scaling_ratios.resizeArray(d_number_levels);
   }

   if (d_is_multiblock) {
      tbox::Pointer< mblk::MultiblockPatchHierarchy<DIM> > mblk_hierarchy =
         hierarchy;
      int nblocks = mblk_hierarchy->getNumberBlocks();
      int level_counter = 0;
      for (int ln = 0; ln <= mblk_hierarchy->getFinestLevelNumber(); ln++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > level =
            mblk_hierarchy->getPatchLevel(ln);
         for (int b = 0; b < nblocks; b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               if (ln == 0) {
                  d_scaling_ratios[level_counter] = hier::IntVector<DIM>(1);
               } else {
                  d_scaling_ratios[level_counter] =
                     patch_level->getRatioToCoarserLevel();
               }
               level_counter++;
            }
         }
      }
   } else {
      for (int ln = 1; ln <= hierarchy->getFinestLevelNumber(); ln++) {
         tbox::Pointer<hier::PatchLevel<DIM> > level =
            hierarchy->getPatchLevel(ln);
         d_scaling_ratios[ln] = level->getRatioToCoarserLevel();
      }
   }

   if (d_top_level_directory_name.empty()) {
      TBOX_ERROR("VisItDataWriter<DIM>::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n     Dump Directory Name is not set" << endl);
   }

   int num_items_to_plot = d_number_visit_variables_plus_depth 
      + d_materials_names.getSize() // number of materials
      + d_number_species;
   
   if (num_items_to_plot == 0) {
      TBOX_ERROR("VisItDataWriter<DIM>::writePlotData"
                 << "\n    No VisIt variables have been registered."
                 << endl);
   }
   
   initializePlotVariableMinMaxInfo(hierarchy);
   
   writeHDFFiles(hierarchy, simulation_time);

}


/*
*************************************************************************
*                                                                       *
* Write plot data from given hierarchy to HDF file                      * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::setSummaryFilename(
   string& filename)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!filename.empty());
#endif
   d_summary_filename = filename + ".samrai";
}

/*
*************************************************************************
*                                                                       *
* Private function to initialize min/max information for the plot       *
* components.  This method will allocate space for the d_mm array on    *
* the VISIT_MASTER processsor (which holds min/max info for all plot    *
* variables on all patches) and will allocate the d_worker_min_max array   *
* on all processors except the VISIT_MASTER.  This latter array is used *
* to store data to be sent to the master when summary information is    *
* written.                                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::initializePlotVariableMinMaxInfo(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif
   
   /*
    * Compute max number of patches on this processor.
    */
   int number_local_patches = 0;
   int tot_number_of_patches = 0;
 
   if (d_is_multiblock) {
      for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
            hierarchy->getPatchLevel(ln);
         for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               tot_number_of_patches += patch_level->getNumberOfPatches();
               for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
                    ip; ip++) {
                  number_local_patches++;
               }
            }
         }
      }
   } else {
      for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
            hierarchy->getPatchLevel(ln);
         tot_number_of_patches += patch_level->getNumberOfPatches();
         for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
              ip; ip++) {
            number_local_patches++;
         }
      }
   }

   int max_number_local_patches = tbox::MPI::maxReduction(number_local_patches);

   /*
    * Determine number of worker processors.  NOTE: subract one because we 
    * don't want to count processor zero.
    */
   int count_me_in = 0;
   if (number_local_patches > 0) count_me_in = 1;
   d_number_working_slaves = tbox::MPI::sumReduction(count_me_in);
   d_number_working_slaves -= 1;
   
   if (tbox::MPI::getRank() != VISIT_MASTER) {
      
      /*
       * Worker processor:  Allocate an array large enough to hold patch
       * min max information.  Pack array by var_item, component number, 
       * level, and local patch number.
       */
      int num_items_to_plot = d_number_visit_variables_plus_depth
         + d_materials_names.getSize() // number materials 
         + d_number_species;

      int num_components = max_number_local_patches*num_items_to_plot;
      if (d_worker_min_max != (patchMinMaxStruct*)NULL) {
         delete [] d_worker_min_max;
      }
      d_worker_min_max = new patchMinMaxStruct[num_components];
      for (int i = 0; i < num_components; i++) {
         d_worker_min_max[i].patch_data_on_disk = false;
         d_worker_min_max[i].min = tbox::IEEE::getDBL_MAX();
         d_worker_min_max[i].max = tbox::IEEE::getDBL_MIN();
         d_worker_min_max[i].material_composition_code = 
            VisMaterialsDataStrategy<DIM>::VISIT_MIXED;
         d_worker_min_max[i].species_composition_code = 
            VisMaterialsDataStrategy<DIM>::VISIT_MIXED;
      }
      
   } else {   // (tbox::MPI::getRank() == VISIT_MASTER)

      /*
       * Master processor:  allocate array for each plot item to hold 
       * min/max information for ALL patches, on all levels. 
       */
      for (typename tbox::List<VisItItem>::Iterator 
              ipi(d_plot_items); ipi; ipi++) {

         for (int comp = 0; comp < ipi().d_depth; comp++) {

            /*
             * Create space for master min/max struct, if it doesn't
             * already exist.
             */
            if (ipi().d_master_min_max[comp] != (patchMinMaxStruct*)NULL) {
               delete [] ipi().d_master_min_max[comp];
            }
            patchMinMaxStruct *mm = 
               new patchMinMaxStruct[tot_number_of_patches];
            ipi().d_master_min_max[comp] = mm;

            for (int pn = 0; pn < number_local_patches; pn++) {
               ipi().d_master_min_max[comp][pn].patch_data_on_disk = false;
               ipi().d_master_min_max[comp][pn].min = tbox::IEEE::getDBL_MAX();
               ipi().d_master_min_max[comp][pn].max = tbox::IEEE::getDBL_MIN();
               ipi().d_master_min_max[comp][pn].material_composition_code = 
                  VisMaterialsDataStrategy<DIM>::VISIT_MIXED;
               ipi().d_master_min_max[comp][pn].species_composition_code = 
                  VisMaterialsDataStrategy<DIM>::VISIT_MIXED;
            }
         }
      }
   } // proc == VISIT_MASTER
   
}

/*
*************************************************************************
*                                                                       *
* Private function to coordinate writing HDF plot files.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::writeHDFFiles(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   double simulation_time)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif

   char temp_buf[VISIT_NAME_BUFSIZE];
   string dump_dirname;
   tbox::HDFDatabase* visit_HDFFilePointer;

   int num_procs = tbox::MPI::getNodes();
   int my_proc = tbox::MPI::getRank();


   if (d_file_cluster_size > num_procs) {
      d_file_cluster_size = num_procs;
   }
   d_my_file_cluster_number = my_proc / d_file_cluster_size;
   d_my_rank_in_file_cluster = my_proc % d_file_cluster_size;

   if (d_my_rank_in_file_cluster == 0) {
      d_file_cluster_leader = true;
   } else {
      d_file_cluster_leader = false;
   }    
   d_number_file_clusters = (int)ceil((double) num_procs / 
                                      (double) d_file_cluster_size);
   d_number_files_this_file_cluster = d_file_cluster_size;

   if (d_my_file_cluster_number == (d_number_file_clusters - 1)) {
      // set d_number_files_this_file_cluster for last cluster
      d_number_files_this_file_cluster = 
        num_procs - ((int) d_file_cluster_size *
           ((int)floor((double) num_procs 
                     / (double) d_file_cluster_size)));
      if (d_number_files_this_file_cluster == 0) {
         d_number_files_this_file_cluster = d_file_cluster_size;
      }
   }

   d_processor_in_file_cluster_number.resizeArray(num_procs);
   for (int i = 0; i < num_procs; i++) {
      d_processor_in_file_cluster_number[i] = i / d_file_cluster_size;
   }

   tbox::Pointer<tbox::HDFDatabase> processor_HDFGroup;
   sprintf(temp_buf, "%05d",d_time_step_number);
   d_current_dump_directory_name = "visit_dump.";
   d_current_dump_directory_name += temp_buf;
   dump_dirname = d_top_level_directory_name + "/";
   dump_dirname = dump_dirname + d_current_dump_directory_name;
   tbox::Utilities::recursiveMkdir(dump_dirname);

   dumpWriteBarrierBegin();
   { 
      // cluster_leader guaranteed to enter this section before anyone else
      sprintf(temp_buf, "/processor_cluster.%05d.samrai", 
              d_my_file_cluster_number);
      string database_name(temp_buf);
      string visit_HDFFilename = dump_dirname + database_name;
      visit_HDFFilePointer = new tbox::HDFDatabase(database_name);
      if (d_file_cluster_leader) {
         visit_HDFFilePointer->mount(visit_HDFFilename, "W"); 
         // "W" creates the HDF file: 
         //      dirname/visit_dump.000n/processor_cluster.000m.samrai 
         //      where n is timestep #, m is processor number
      } else { // file already created other procs just need to open it
         if (visit_HDFFilePointer->mount(visit_HDFFilename, "R") < 0){
            TBOX_ERROR("VisItDataWriter<DIM>::writeHDFFiles"
               << "\n    data writer with name " << d_object_name
               << "\n    Error attempting to open visit file " 
               <<        visit_HDFFilename << endl);
         }
      }

      // create group for this proc
      sprintf(temp_buf, "processor.%05d", my_proc);
      processor_HDFGroup = 
           visit_HDFFilePointer->putDatabase(string(temp_buf));
      writeVisItVariablesToHDFFile(processor_HDFGroup,
                                   hierarchy,
                                   0, 
                                   hierarchy->getFinestLevelNumber());
      visit_HDFFilePointer->unmount(); // invokes H5FClose
      delete visit_HDFFilePointer; // deletes tbox::HDFDatabase object
   }   

   dumpWriteBarrierEnd();

   tbox::MPI::barrier();

   writeSummaryToHDFFile(dump_dirname, 
                         hierarchy, 
                         0, 
                         hierarchy->getFinestLevelNumber(), 
                         simulation_time);
}

/*
*************************************************************************
*                                                                       *
* Private function to find global patch number, given level number and  *
* local patch number. This is needed because SAMRAI maintains the patch *
* number on each level, but VisIt needs a unique number for each patch  *
* written.                                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VisItDataWriter<DIM>::getGlobalPatchNumber(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const int block_number,
   const int patch_number)
{
   int global_patch_id=0;

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(level_number >= 0);
   assert(block_number >= -1);
   assert(patch_number >= 0);
#endif

   if (d_is_multiblock) {
      for (int i = 0; i < level_number; i++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
            hierarchy->getPatchLevel(i);
         int nblocks = mblk_level->getNumberBlocks();
         for (int b = 0; b < nblocks; b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               global_patch_id += patch_level->getNumberOfPatches();
            }
         }
      }
      tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > this_mblk_level =
         hierarchy->getPatchLevel(level_number);

      for (int b = 0; b < block_number; b++) {
         tbox::Pointer< hier::PatchLevel<DIM> > this_patch_level =
            this_mblk_level->getPatchLevelForBlock(b);
         if (!this_patch_level.isNull()) {
            global_patch_id += this_patch_level->getNumberOfPatches();
         }
      }

   } else {

      for (int i = 0; i < level_number; i++) {
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
            hierarchy->getPatchLevel(i);
         global_patch_id += patch_level->getNumberOfPatches();
      }

   }

   global_patch_id += patch_number;
   return(global_patch_id);
}

/*
*************************************************************************
*                                                                       *
* Private function to write variables & materials data to an HDF File.  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::writeVisItVariablesToHDFFile(
   tbox::Pointer<tbox::HDFDatabase> processor_HDFGroup,
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   int coarsest_level,
   int finest_level)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(coarsest_level >= 0);
   assert(finest_level >= 0);
#endif

   /*
    * Reset the var_id_ctr - this is used to record min/max summary 
    * information for every plotted variable on the patch.  It is incremented
    * for each component (i.e. depth), of each variable, of each patch.
    */
   d_var_id_ctr = 0;
   
   char temp_buf[VISIT_NAME_BUFSIZE];
   tbox::Pointer<tbox::HDFDatabase> level_HDFGroup, patch_HDFGroup;

   if (d_is_multiblock) {

      int unique_level_number = coarsest_level;
      for (int ln = coarsest_level; ln <= finest_level; ln++) {

         /*
          * create new HDFGroup for this level
          */
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
                                            hierarchy->getPatchLevel(ln);

         int nblocks = mblk_level->getNumberBlocks();

         int patch_counter = 0;
         for (int b = 0; b < nblocks; b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);

            if (!patch_level.isNull()) {

               sprintf(temp_buf, "level.%05d", unique_level_number);
               level_HDFGroup =
                  processor_HDFGroup->putDatabase(string(temp_buf));
               unique_level_number++;

               hier::IntVector<DIM> coarsen_ratio =
                  patch_level->getRatioToCoarserLevel();

               for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
                    ip; ip++) {
                  tbox::Pointer< hier::Patch<DIM> > patch =
                     patch_level->getPatch(ip());

                  /*
                   * create new HDFGroup for this patch
                   */
                  int pn = patch->getPatchNumber();
                  sprintf(temp_buf, "patch.%05d",pn);
                  patch_HDFGroup =
                     level_HDFGroup->putDatabase(string(temp_buf));

                  packRegularAndDerivedData(patch_HDFGroup,
                                            hierarchy,
                                            ln,
                                            b,
                                            *patch);

                  if (d_materials_names.getSize() > 0) {
                     packMaterialsData(patch_HDFGroup,
                                       hierarchy,
                                       ln,
                                       b,
                                       *patch);

                     packSpeciesData(hierarchy,
                                     ln,
                                     b,
                                     *patch);
                  }

               }
               patch_counter += patch_level->getNumberOfPatches();
            }
         }
      }
   } else {
      for (int ln = coarsest_level; ln <= finest_level; ln++) {

         /*
          * create new HDFGroup for this level
          */
         sprintf(temp_buf, "level.%05d", ln);
         level_HDFGroup = processor_HDFGroup->putDatabase(string(temp_buf));

         tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
                                         hierarchy->getPatchLevel(ln);
         hier::IntVector<DIM> coarsen_ratio =
              patch_level->getRatioToCoarserLevel();

         for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
              ip; ip++) {
            tbox::Pointer< hier::Patch<DIM> > patch =
               patch_level->getPatch(ip());

            /*
             * create new HDFGroup for this patch
             */
            int pn = patch->getPatchNumber();
            sprintf(temp_buf, "patch.%05d",pn);
            patch_HDFGroup = level_HDFGroup->putDatabase(string(temp_buf));

            int bn = -1;
            packRegularAndDerivedData(patch_HDFGroup,
                                      hierarchy,
                                      ln,
                                      bn,
                                      *patch);
         
            if (d_materials_names.getSize() > 0) {
               packMaterialsData(patch_HDFGroup,
                                 hierarchy,
                                 ln,
                                 bn,
                                 *patch);

               packSpeciesData(hierarchy,
                               ln,
                               bn,
                               *patch);
            }

         }
      }
   }

   /* 
    * Clean up from packing operations.  If this is not done there is a dangling smart
    * pointer reference to HDF5 groups and the file may not be written/closed.
    */
   for (typename tbox::List<VisItItem>::Iterator 
	   ipi(d_plot_items); ipi; ipi++) {
      ipi().d_species_HDFGroup = NULL;
      ipi().d_extents_species_HDFGroup = NULL;
   }
}

/*
*************************************************************************
*                                                                       *
* Private function to pack regular and derived VisIt variables into     *
* specified HDF database.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::packRegularAndDerivedData(
   tbox::Pointer<tbox::HDFDatabase> patch_HDFGroup,
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const int block_number,
   hier::Patch<DIM>& patch)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(level_number >= 0);
   assert(block_number >= -1);
#endif
   
   /*
    * Loop over variables and write out those that are NOT
    * material or species variables.
    */
   for (typename tbox::List<VisItItem>::Iterator 
           ipi(d_plot_items); ipi; ipi++) {

      /*
       * Only write regular (non-derived) and derived vars, not
       * material or spacies vars.
       */
      if (!(ipi().d_isa_material || ipi().d_isa_species)) {
      
         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(), 
                                      hier::IntVector<DIM>(0),
                                      ipi().d_var_centering);
         
         double *dbuffer = new double[buf_size]; // used to pack var
         float *fbuffer = new float[buf_size]; // copy to float for writing

         for (int depth_id = 0; depth_id < ipi().d_depth; depth_id++) {
            
            /*
             * If its derived data, pack via the derived writer.
             * Otherwise, pack with local private method.
             */
            bool data_exists_on_patch = false;
            int patch_data_id = VISIT_UNDEFINED_INDEX;
            if (ipi().d_is_derived) {
               
               // derived data
               data_exists_on_patch = 
                  ipi().d_derived_writer->
                  packDerivedDataIntoDoubleBuffer(
                     dbuffer,
                     patch,
                     patch.getBox(),
                     ipi().d_var_name,
                     depth_id);

            } else {
               
               /*
                * Check if patch data id has been reset on the level.  If
                * not, just use the original registered data id.
                */
               patch_data_id = ipi().d_patch_data_index;
               if (ipi().d_level_patch_data_index.getSize() > level_number) {
                  patch_data_id =
                     ipi().d_level_patch_data_index[level_number];
               }
               
               data_exists_on_patch =
                  patch.checkAllocated(patch_data_id);
               
               if (data_exists_on_patch) {
                  int new_depth_id = ipi().d_start_depth_index + depth_id;
                  
                  // regular (non-derived) data
                  packPatchDataIntoDoubleBuffer(
                     patch.getPatchData(patch_data_id),
                     new_depth_id,
                     ipi().d_var_data_type,
                     patch.getBox(),
                     dbuffer,
                     ipi().d_var_centering);
               }
            }                  

            double dmax = tbox::IEEE::getDBL_MIN();
            double dmin = tbox::IEEE::getDBL_MAX();

            if (data_exists_on_patch) {

               /*
                * Scale data (while still double)
                */
               const double scale = ipi().d_scale_factor;
               if (scale != 1.0) {
                  for (int i = 0; i < buf_size; i++) {
                     dbuffer[i] *= scale;
                  }
               }

               /*
                * Determine patch min/max.
                */

               int i;
               for (i = 0; i < buf_size; i++) {
                  if (dbuffer[i] > dmax) {
                     dmax = dbuffer[i];
                  }
                  if (dbuffer[i] < dmin) {
                     dmin = dbuffer[i];
                  }
               }

               checkFloatMinMax(dmin,
                                dmax,
                                level_number,
                                patch.getPatchNumber(),
                                patch_data_id);
               
               /*
                * Convert buffer from double to float
                */
               for (i = 0; i < buf_size; i++) {
                  fbuffer[i] = (float) dbuffer[i];
               }
               
               /*
                * Write to disk
                */
               string vname = ipi().d_visit_var_name[depth_id];
               patch_HDFGroup->putFloatArray(vname,
                                             fbuffer,
                                             buf_size);
                  
            } else { // data does not exist on patch

               dmax = 0.;
               dmin = 0.;

            }
            
             
            /*
             * Write min/max summary info
             */
            if (tbox::MPI::getRank() == VISIT_MASTER) {
               int gpn = getGlobalPatchNumber(hierarchy,
                                              level_number,
                                              block_number,
                                              patch.getPatchNumber());
               ipi().d_master_min_max[depth_id][gpn].patch_data_on_disk = 
                  data_exists_on_patch;
               ipi().d_master_min_max[depth_id][gpn].min = dmin;
               ipi().d_master_min_max[depth_id][gpn].max = dmax;
            } else {
               d_worker_min_max[d_var_id_ctr].patch_data_on_disk = 
                  data_exists_on_patch;
               d_worker_min_max[d_var_id_ctr].min = dmin;
               d_worker_min_max[d_var_id_ctr].max = dmax;
            }
            

            /*
             * Increment local var_id counter used for d_mm array.
             */
            d_var_id_ctr++;
            
         } // loop over var depths

         delete [] dbuffer;
         delete [] fbuffer;
         
      } // var is not species or material
      
   } // iterate over vars
   
}


/*
*************************************************************************
*                                                                       *
* Private function to pack Material variables into specified HDF        *
* database.                                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::packMaterialsData(
   tbox::Pointer<tbox::HDFDatabase> patch_HDFGroup,
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const int block_number,
   hier::Patch<DIM>& patch)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(level_number >= 0);
   assert(block_number >= -1);
#endif


   
   /*
    * Loop over variables and pull out those that are material variables.
    */
   tbox::Pointer<tbox::HDFDatabase> materials_HDFGroup;
   tbox::Pointer<tbox::HDFDatabase> material_name_HDFGroup;
   for (typename tbox::List<VisItItem>::Iterator 
           ipi(d_plot_items); ipi; ipi++) {

      if (ipi().d_isa_material) {
      
         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(), 
                                      hier::IntVector<DIM>(0),
                                      ipi().d_var_centering);
         
         double *dbuffer = new double[buf_size]; // used to pack var
         float *fbuffer = new float[buf_size]; // copy to float for writing

         for (int depth_id = 0; depth_id < ipi().d_depth; depth_id++) {
            
            /*
             * Create the "materials" HDF database entry:
             * <patch HDF group>
             *    "materials"  (only if there are materials)
             *       material_name
             *          "species"    (only if matl has species)
             */
            if (materials_HDFGroup.isNull()) {
               materials_HDFGroup = 
                  patch_HDFGroup->putDatabase("materials");
            }  
            
            // create material_name HDF database
            string mname = ipi().d_material_name;
            material_name_HDFGroup = 
               materials_HDFGroup->putDatabase(mname);
            
            // create "species" HDF database for material name
            if (ipi().d_species_names.getSize() > 0) {
               ipi().d_species_HDFGroup = 
                  material_name_HDFGroup -> putDatabase("species");
            }
            
            // pack the buffer with material data
            int return_code = d_materials_writer->
               packMaterialFractionsIntoDoubleBuffer(
                  dbuffer,
                  patch,
                  patch.getBox(),
                  ipi().d_material_name);
         
            
            // check return code
            if ((return_code != VisMaterialsDataStrategy<DIM>::VISIT_MIXED) 
                && (return_code != 
                    VisMaterialsDataStrategy<DIM>::VISIT_ALLONE)
                && (return_code != 
                    VisMaterialsDataStrategy<DIM>::VISIT_ALLZERO)) { 
               TBOX_ERROR(
                  "VisItDataWriter<DIM>::packMaterialsData()"
                  << "\n    Invalid return value from "
                  << "packMaterialFractionsIntoDoubleBuffer()\n");
            }            

            double dmax = tbox::IEEE::getDBL_MIN();
            double dmin = tbox::IEEE::getDBL_MAX();
            bool data_on_disk = false;
            
            if (return_code == VisMaterialsDataStrategy<DIM>::VISIT_MIXED) {

               data_on_disk = true;
               
               /*
                * Scale data (while still double)
                */
               const double scale = ipi().d_scale_factor;
               if (scale != 1.0) {
                  for (int i = 0; i < buf_size; i++) {
                     dbuffer[i] *= scale;
                  }
               }

               /*
                * Determine patch min/max.
                */
               int i;
               for (i = 0; i < buf_size; i++) {
                  if (dbuffer[i] > dmax) {
                     dmax = dbuffer[i];
                  }
                  if (dbuffer[i] < dmin) {
                     dmin = dbuffer[i];
                  }
               }

               int dummy_pdata_id = VISIT_UNDEFINED_INDEX;
               checkFloatMinMax(dmin,
                                dmax,
                                level_number,
                                patch.getPatchNumber(),
                                dummy_pdata_id);

               /*
                * Convert buffer from double to float
                */
               for (i = 0; i < buf_size; i++) {
                  fbuffer[i] = (float) dbuffer[i];
               }
               
               /*
                * Write to disk
                */
               string vname = ipi().d_material_name;
               vname = vname + "-fractions";
               material_name_HDFGroup -> putFloatArray(vname,
						       fbuffer,
						       buf_size);
                  
            } else if (return_code == 
                       VisMaterialsDataStrategy<DIM>::VISIT_ALLONE) {
 
               data_on_disk = false;
               dmin = 1.0;
               dmax = 1.0;

            } else { // return code == VISIT_ALLZERO  

               data_on_disk = false;
               dmin = 0.0;
               dmax = 0.0;

            }

            /*
             * Write min/max summary info
             */
            if (tbox::MPI::getRank() == VISIT_MASTER) {
               int gpn = getGlobalPatchNumber(hierarchy,
                                              level_number,
                                              block_number,
                                              patch.getPatchNumber());
               ipi().d_master_min_max[depth_id][gpn].patch_data_on_disk = 
                  data_on_disk;
               ipi().d_master_min_max[depth_id][gpn].min = dmin;
               ipi().d_master_min_max[depth_id][gpn].max = dmax;
               ipi().d_master_min_max[depth_id][gpn].
                  material_composition_code = return_code;
            } else {
               d_worker_min_max[d_var_id_ctr].patch_data_on_disk = 
                  data_on_disk;
               d_worker_min_max[d_var_id_ctr].min = dmin;
               d_worker_min_max[d_var_id_ctr].max = dmax;
               d_worker_min_max[d_var_id_ctr].material_composition_code = 
                  return_code;
            }

            /*
             * Increment local var_id counter used for d_mm array.
             */
            d_var_id_ctr++;
            
         } // loop over var depths
         
         delete [] fbuffer;
         delete [] dbuffer;
      
      } // var is a material

   } // iterate over vars
   
}

/*
*************************************************************************
*                                                                       *
* Private function to pack Species variables into specified HDF         *
* database.                                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::packSpeciesData(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const int block_number,
   hier::Patch<DIM>& patch)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(level_number>=0);
   assert(block_number>=-1);   
#endif

   /*
    * Loop over variables and pull out those that are material variables.
    */
   for (typename tbox::List<VisItItem>::Iterator 
           ipi(d_plot_items); ipi; ipi++) {

      if (ipi().d_isa_species) {
      
         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(), 
                                      hier::IntVector<DIM>(0),
                                      ipi().d_var_centering);
         
         double *dbuffer = new double[buf_size]; // used to pack var
         float *fbuffer = new float[buf_size]; // copy to float for writing

         for (int depth_id = 0; depth_id < ipi().d_depth; depth_id++) {
            
            // pack the buffer with species data
            int return_code = d_materials_writer->
               packSpeciesFractionsIntoDoubleBuffer(
                  dbuffer,
                  patch,
                  patch.getBox(),
                  ipi().d_material_name,
                  ipi().d_species_name);
         
            
            // check return code
            if ((return_code != VisMaterialsDataStrategy<DIM>::VISIT_MIXED) 
                && (return_code != 
                    VisMaterialsDataStrategy<DIM>::VISIT_ALLONE)
                && (return_code == 
                    VisMaterialsDataStrategy<DIM>::VISIT_ALLZERO)) { 
               TBOX_ERROR("VisItDataWriter<DIM>::packSpeciesData()"
                          << "\n    Invalid return value from "
                          << "packSpeciesFractionsIntoDoubleBuffer()\n"
                          << endl);
            }            

            double dmax = tbox::IEEE::getDBL_MIN();
            double dmin = tbox::IEEE::getDBL_MAX();
            bool data_on_disk = false;

            if (return_code == VisMaterialsDataStrategy<DIM>::VISIT_MIXED) {

               /*
                * Scale data (while still double)
                */
               const double scale = ipi().d_scale_factor;
               if (scale != 1.0) {
                  for (int i = 0; i < buf_size; i++) {
                     dbuffer[i] *= scale;
                  }
               }

               /*
                * Determine patch min/max.
                */
               int i;
               for (i = 0; i < buf_size; i++) {
                  if (dbuffer[i] > dmax) {
                     dmax = dbuffer[i];
                  }
                  if (dbuffer[i] < dmin) {
                     dmin = dbuffer[i];
                  }
               }

               int dummy_pdata_id = VISIT_UNDEFINED_INDEX;
               checkFloatMinMax(dmin,
                                dmax,
                                level_number,
                                patch.getPatchNumber(),
                                dummy_pdata_id);

               /*
                * Convert buffer from double to float
                */
               for (i = 0; i < buf_size; i++) {
                  fbuffer[i] = (float) dbuffer[i];
               }
               
               /*
                * Write to disk
                */
               string sname = ipi().d_species_name;
               ipi().d_parent_material_pointer->
                  d_species_HDFGroup->putFloatArray(sname,
                                                    fbuffer,
                                                    buf_size);

            } else if (return_code == 
                       VisMaterialsDataStrategy<DIM>::VISIT_ALLONE) {
 
               data_on_disk = false;
               dmin = 1.0;
               dmax = 1.0;

            } else { // return code == VISIT_ALLZERO  

               data_on_disk = false;
               dmin = 0.0;
               dmax = 0.0;

            }

            /*
             * Write min/max summary info
             */
            if (tbox::MPI::getRank() == VISIT_MASTER) {
               int gpn = getGlobalPatchNumber(hierarchy,
                                              level_number,
                                              block_number,
                                              patch.getPatchNumber());
               ipi().d_master_min_max[depth_id][gpn].patch_data_on_disk = 
                  data_on_disk;
               ipi().d_master_min_max[depth_id][gpn].min = dmin;
               ipi().d_master_min_max[depth_id][gpn].max = dmax;
               ipi().d_master_min_max[depth_id][gpn].species_composition_code = 
                  return_code;
            } else {
               d_worker_min_max[d_var_id_ctr].patch_data_on_disk = 
                  data_on_disk;
               d_worker_min_max[d_var_id_ctr].min = dmin;
               d_worker_min_max[d_var_id_ctr].max = dmax;
               d_worker_min_max[d_var_id_ctr].species_composition_code = 
                  return_code;
            }

            /*
             * Increment local var_id counter used for d_mm array.
             */
            d_var_id_ctr++;
            
         } // loop over var depths

         delete [] fbuffer;
         delete [] dbuffer;
         
      } // var is a species
      
   } // iterate over vars
   
}


/*
*************************************************************************
*                                                                       *
* Private function to check float min/max values.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::checkFloatMinMax(
   const double dmin,
   const double dmax,
   const int level_number,
   const int patch_number,
   const int patch_data_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number>=0);
   assert(patch_number>=0);
   assert(patch_data_id>=-1);   
#endif

   double fmin = -(double)tbox::IEEE::getFLT_MAX();
   double fmax = (double)tbox::IEEE::getFLT_MAX();
   
   if (dmin < fmin) {
      TBOX_ERROR("VisItDataWriter<DIM>:"
                 << "\n    hier::Patch data is less than FLT_MIN "
                 << "\n    level: " << level_number 
                 <<"  patch: " << patch_number
                 <<"  patch_data_id: " << patch_data_id
                 <<"  value: " << dmin
                 <<"\n    It cannot be read by VisIt." 
                 <<"\n    Make sure data is properly initialized or" 
                 <<"\n    use scale factor to increase its size.");
   }
   if (dmax > fmax) {
      TBOX_ERROR("VisItDataWriter<DIM>:"
                 << "\n    hier::Patch data is greater than FLT_MAX "
                 << "\n    level: " << level_number 
                 <<"  patch: " << patch_number
                 <<"  patch_data_id: " << patch_data_id
                 <<"  value: " << dmax
                 <<"\n    It cannot be interpreted by VisIt." 
                 <<"\n    Make sure data is properly initialized or" 
                 <<"\n    use scale factor to decrease its size.");
   }            
}

/*
*************************************************************************
*                                                                       *
* Private function to write one summary HDF file covering data from all *
* processors for use by VisIt.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::writeSummaryToHDFFile(
   string dump_dirname,
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   int coarsest_plot_level,
   int finest_plot_level,
   double simulation_time)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(coarsest_plot_level>=0);
   assert(finest_plot_level>=0);
#endif

   int i, ln, pn;
   
   /*
    * Pack patch min/max information 
    */
   exchangeMinMaxPatchInformation(hierarchy, 
                                  coarsest_plot_level, 
                                  finest_plot_level);

   /*
    * The "VISIT_MASTER" writes a set of summary information to 
    * the summary file that describes data contained in the visit 
    * files written by each tbox::MPI process. 
    */
   int my_proc = tbox::MPI::getRank();
   if (my_proc == VISIT_MASTER) {
      char temp_buf[VISIT_NAME_BUFSIZE];
      //sprintf(temp_buf, "/summary.samrai");
      //string summary_HDFFilename = dump_dirname + temp_buf;
      string summary_HDFFilename = dump_dirname + "/" + d_summary_filename;
      tbox::Pointer<tbox::HDFDatabase> summary_HDFFilePointer = 
         new tbox::HDFDatabase("root");
      summary_HDFFilePointer->mount(summary_HDFFilename, "W");

      /*
       * Create BASIC information HDF Group and provide it the following 
       * information:
       *   - VisItWriter version number
       *   - grid_type (CARTESIAN or DEFORMED)
       *   - simulation time
       *   - time step number
       *   - number of processors
       *   - number of file clusters
       *   - DIM
       *   - number of levels
       *   - number of patches at each level (array - int[nlevels])
       *   - total number of patches (on all levels)
       *   - ratio to coarser level (array - int[nlevels][ndim])
       */

      sprintf(temp_buf, "BASIC_INFO");
      tbox::Pointer<tbox::HDFDatabase> basic_HDFGroup
                  = summary_HDFFilePointer->putDatabase(string(temp_buf));
      hid_t basic_group_id = basic_HDFGroup->getGroupId();

      string key_string = "VDR_version_number";
      basic_HDFGroup->putFloat(key_string, VISIT_DATAWRITER_VERSION_NUMBER);

      key_string = "grid_type";
      string data_string;
      if (d_grid_type == VISIT_CARTESIAN) {
         data_string = "CARTESIAN";
      } else if (d_grid_type == VISIT_DEFORMED) {
          data_string = "DEFORMED";
      } else {
         TBOX_ERROR("VisItDataWriter<DIM>::writeSummaryToHDFFile"
            << "\n    data writer with name " << d_object_name
            << "\n    Illegal grid type: " << d_grid_type << endl);
      }
      basic_HDFGroup->putString(key_string, data_string);

      key_string = "time";
      basic_HDFGroup->putDouble(key_string, simulation_time);

      key_string = "time_step_number";
      basic_HDFGroup->putInteger(key_string, d_time_step_number);

      key_string = "number_processors";
      basic_HDFGroup->putInteger(key_string, tbox::MPI::getNodes());

      key_string = "number_file_clusters";
      basic_HDFGroup->putInteger(key_string, d_number_file_clusters);

      key_string = "number_dimensions_of_problem";
      basic_HDFGroup->putInteger(key_string, DIM);

      int num_levels;
      int tot_number_of_patches = 0;
      if (d_is_multiblock) {
         tbox::Pointer< mblk::MultiblockPatchHierarchy<DIM> > mblk_hierarchy =
            hierarchy;

         //key_string = "number_levels";
         num_levels =
            (mblk_hierarchy->getNumberLevels()) *
            (mblk_hierarchy->getNumberBlocks());

         //key_string = "number_patches_at_level";
         tbox::Array<int> num_patches_per_level(num_levels);
         int unique_level_number = 0;
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            num_patches_per_level[ln] = 0;
            tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
               hierarchy->getPatchLevel(ln);
            for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
               tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
                  mblk_level->getPatchLevelForBlock(b);
               if (!patch_level.isNull()) {
                  num_patches_per_level[unique_level_number] =
                     patch_level->getNumberOfPatches();
                  tot_number_of_patches += patch_level->getNumberOfPatches();
                  unique_level_number++;
               }
            }
         }

         num_levels = unique_level_number;
         num_patches_per_level.resizeArray(num_levels);

         key_string = "number_levels";
         basic_HDFGroup->putInteger(key_string, num_levels);

         key_string = "number_patches_at_level";
         basic_HDFGroup->putIntegerArray(key_string,
                                         num_patches_per_level);

         key_string = "number_global_patches";
         basic_HDFGroup->putInteger(key_string, tot_number_of_patches);

      } else {

         key_string = "number_levels";
         num_levels = hierarchy->getNumberLevels();
         basic_HDFGroup->putInteger(key_string, num_levels);

         key_string = "number_patches_at_level";
         tbox::Array<int> num_patches_per_level(num_levels);
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
               hierarchy->getPatchLevel(ln);
            num_patches_per_level[ln] = patch_level->getNumberOfPatches();
            tot_number_of_patches += num_patches_per_level[ln];
         }
         basic_HDFGroup->putIntegerArray(key_string,
                                         num_patches_per_level);

         key_string = "number_global_patches";
         basic_HDFGroup->putInteger(key_string, tot_number_of_patches);
      }

      /*
       * When writing Visit data, it expects to see 3D data for 
       * xlo, dx, ratios_to_coarser, and number of ghosts.  The 
       * VISIT_FIXED_DIM is set to 3, and the third element
       * is zero if we have 2D data.
       */
      key_string = "ratios_to_coarser_levels";
      int idx = 0;
      int* rtcl = new int[num_levels * VISIT_FIXED_DIM];
      for (i = 0; i < num_levels * VISIT_FIXED_DIM; i++) rtcl[i] = 0;
      for (ln = 0; ln <= finest_plot_level; ln++) {
         for (i = 0; i < DIM; i++) { 
            idx = ln * VISIT_FIXED_DIM + i;
            if (ln == 0) {
               rtcl[idx] = VISIT_UNDEFINED_INDEX;
            } else {
               rtcl[idx] = d_scaling_ratios[ln](i);
            }
         }
      }
      HDFputIntegerArray2D(key_string,
                           rtcl, 
                           num_levels, 
                           VISIT_FIXED_DIM,
                           basic_group_id);
      delete [] rtcl;

      /*
       * Write to BASIC HDF group information about the
       * VisIt plot variables:
       *   - number of visit variables
       *
       *   for each variable {
       *      - name
       *      - centering (1 = CELL, 0 = NODE)
       *      - scale factor
       *      - depth
       *      - ghosts
       *   }
       */
      key_string = "number_visit_variables";
      basic_HDFGroup->putInteger(key_string, d_number_visit_variables);

      tbox::Array<string> var_names(d_number_visit_variables);
      tbox::Array<int> var_centering(d_number_visit_variables);
      tbox::Array<double> var_scale_factors(d_number_visit_variables);
      tbox::Array<int> var_depths(d_number_visit_variables);
      int* var_ghosts = new int[d_number_visit_variables * VISIT_FIXED_DIM];
      for (i = 0; i < d_number_visit_variables*VISIT_FIXED_DIM; i++) {
         var_ghosts[i] = 0;
      }
      
      i = 0;
      for (typename tbox::List<VisItItem>::Iterator 
              ipi(d_plot_items); ipi; ipi++) {

         if (!(ipi().d_isa_material || ipi().d_isa_species)) {               
            var_names[i] = ipi().d_var_name;
            if ((ipi().d_var_centering == VISIT_CELL) || 
                (ipi().d_var_centering == VISIT_UNKNOWN_CELL)) {
               var_centering[i] = 1;
            } else {
               var_centering[i] = 0;
            }
            if (ipi().d_is_deformed_coords) {
               var_scale_factors[i] = 0.0;
            } else {
               var_scale_factors[i] = ipi().d_scale_factor;
            }
            var_depths[i] = ipi().d_depth;
            for (int dim = 0; dim < DIM; dim++) {
               var_ghosts[i*VISIT_FIXED_DIM+dim] = 
                  0;
            }
            i++;
         }
      }

      key_string = "var_names";
      basic_HDFGroup->putStringArray(key_string,
                                     var_names);

      key_string = "var_cell_centered";
      basic_HDFGroup->putIntegerArray(key_string,
                                      var_centering);

      key_string = "scaling";
      basic_HDFGroup->putDoubleArray(key_string,
                                     var_scale_factors);

      if (d_grid_type == VISIT_DEFORMED) {
         tbox::Array<double> coord_scaling(VISIT_FIXED_DIM);
         for (typename tbox::List<VisItItem>::Iterator 
                 ipi(d_plot_items); ipi; ipi++) {
            if (ipi().d_is_deformed_coords) {
               for (i = 0; i < VISIT_FIXED_DIM; i++) {
                  coord_scaling[i] = 0.0;
                  if (i < DIM) {
                     coord_scaling[i] = ipi().d_coord_scale_factor[i];
                  }
               }
            }
         }
         key_string = "deformed_coordinate_scaling";
         basic_HDFGroup->putDoubleArray(key_string,
                                        coord_scaling);
      }
      


      key_string = "var_number_components";
      basic_HDFGroup->putIntegerArray(key_string,
                                      var_depths);

      key_string = "var_number_ghosts";
      HDFputIntegerArray2D(key_string,
                          var_ghosts,
                          d_number_visit_variables,
                          VISIT_FIXED_DIM,
                          basic_group_id);
      delete [] var_ghosts;

      /*
       * Write time and data info to BASIC HDF group
       */
      key_string = "time_of_dump";
      string temp = dump_dirname + "/temp";
      ofstream dfile(temp.c_str());
      char cmdstr[100];
      sprintf(cmdstr,"date > %s",temp.c_str());
      system(cmdstr);
      ifstream date_file(temp.c_str());
      string ds1, ds2, ds3, ds4, ds5, ds6;
      date_file >> ds1 >> ds2 >> ds3 >> ds4 >> ds5 >> ds6;
      string date = ds1 + " " + ds2 + " " + ds3 + " " + ds4 
                        + " " + ds5 + " " + ds6;
      basic_HDFGroup->putString(key_string, date);
      dfile.close(); 
      date_file.close();
      sprintf(cmdstr,"rm %s",temp.c_str());
      system(cmdstr);



      /*
       * Write general materials and species information to the 
       * materials and species HDF groups:
       *
       *  - material names
       *  - species names
       */
      tbox::Pointer<tbox::HDFDatabase> materials_HDFGroup;
      tbox::Pointer<tbox::HDFDatabase> species_HDFGroup;
      if (d_materials_names.getSize() > 0) {
         sprintf(temp_buf, "materials");
         materials_HDFGroup = 
                  summary_HDFFilePointer->putDatabase(string(temp_buf));

         key_string = "material_names";
         materials_HDFGroup->putStringArray(key_string,
                                            d_materials_names);

         int mat_ghosts[VISIT_FIXED_DIM];
         for (i = 0; i < VISIT_FIXED_DIM; i++) mat_ghosts[i] = 0;
         for (typename tbox::List<VisItItem>::Iterator 
                 ipi(d_plot_items); ipi; ipi++) {            
            if (ipi().d_isa_material) {
               for (int dim = 0; dim < DIM; dim++) {
                  mat_ghosts[dim] = 
                     0;
               }
            }
         }
         key_string = "material_number_ghosts";
         materials_HDFGroup->putIntegerArray(key_string, 
                                             mat_ghosts,
                                             VISIT_FIXED_DIM);

         sprintf(temp_buf, "species");
         species_HDFGroup = 
                  materials_HDFGroup->putDatabase(string(temp_buf));

         for (i = 0; i < d_materials_names.getSize(); i++) {
            key_string = d_materials_names[i];
            for (typename tbox::List<VisItItem>::Iterator 
                                     ipi(d_plot_items); ipi; ipi++) {
               if ((ipi().d_material_name == d_materials_names[i]) &&
                                            ipi().d_isa_material) {
                  if (ipi().d_species_names.getSize() > 0) {
                     species_HDFGroup->putStringArray(
                        key_string,
                        ipi().d_species_names);
                  }
               }
            }
         }
      }
      
      /*
       * When writing Visit data, it expects to see 3D data for 
       * xlo, dx, ratios_to_coarser, and number of ghosts.  The 
       * VISIT_FIXED_DIM is set to 3, and the third element
       * is zero if we have 2D data.
       */
      double geom_lo[VISIT_FIXED_DIM] = {0.,0.,0.}; 

      double dx_curr_lev[DIM];
      double patch_xlo, patch_xhi;

      /*
       * Add mesh dx information to BASIC group
       */

      double *dx = new double[VISIT_FIXED_DIM * num_levels];
      double *xlo = new double[VISIT_FIXED_DIM * num_levels];
      for (i = 0; i < VISIT_FIXED_DIM*num_levels; i++) {
         dx[i] = 0.0;
         xlo[i] = 0.0;
      }
      if (d_grid_type != VISIT_DEFORMED) {
         //This is never entered in multiblock case
         tbox::Pointer< hier::PatchHierarchy<DIM> > patch_hierarchy =
            hierarchy;
         const tbox::Pointer< geom::CartesianGridGeometry<DIM> > ggeom = 
            patch_hierarchy->getGridGeometry();
         int next = 0;
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            for (i = 0; i < VISIT_FIXED_DIM; i++) {
               if (i < DIM) {
                  if (ln == 0) {
                     xlo[i] = ggeom->getXLower()[i];
                     dx_curr_lev[i] = ggeom->getDx()[i]; // coarsest level dx
                     dx[next] = dx_curr_lev[i];
                  } else {
                     double scale_ratio = (double)d_scaling_ratios[ln](i);
                     dx_curr_lev[i] = dx_curr_lev[i]/scale_ratio;
                     dx[next] = dx_curr_lev[i];
                  }
               }
               next++;
            }
         }
      }
      
      key_string = "dx";
      HDFputDoubleArray2D(key_string,
                          dx, 
                          num_levels, 
                          VISIT_FIXED_DIM,
                          basic_group_id);
      delete [] dx;

      /*
       * Add mesh xlo information to BASIC group
       */
      key_string = "XLO";
      basic_HDFGroup->putDoubleArray(key_string, 
                                     xlo, 
                                     VISIT_FIXED_DIM);
      delete [] xlo;

      /*
       * Write parent/child information
       */
      writeParentChildInfoToSummaryHDFFile(hierarchy, basic_HDFGroup);

      /*
       * Write processor mapping information and domain extents
       * for each patch.
       */

      tbox::Pointer<tbox::HDFDatabase> extents_HDFGroup;
      sprintf(temp_buf, "extents");
      extents_HDFGroup = 
         summary_HDFFilePointer->putDatabase(string(temp_buf));
      hid_t extents_group_id = extents_HDFGroup->getGroupId();

      /*
       * Create "patch_map" sub-database of extents group.
       */
      patchMapStruct *pms = new patchMapStruct[tot_number_of_patches];

      if (d_is_multiblock) {
         int unique_level_number = 0;
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
                                   hierarchy->getPatchLevel(ln);
            for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
               tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
                  mblk_level->getPatchLevelForBlock(b);
               if (!patch_level.isNull()) {
                  tbox::Array<int> proc_mapping =
                     patch_level->getProcessorMapping().getProcessorMapping();

                  for (pn = 0; pn < patch_level->getNumberOfPatches(); pn++) {
                     int proc_num = proc_mapping[pn];
                     int global_patch_id =
                        getGlobalPatchNumber(hierarchy,ln,b,pn);
                     pms[global_patch_id].processor_number = proc_num;
                     pms[global_patch_id].file_cluster_number =
                        d_processor_in_file_cluster_number[proc_num];
                     pms[global_patch_id].level_number = unique_level_number;
                     pms[global_patch_id].patch_number = pn;
                  }
                  unique_level_number++;
               }
            }
         }
      } else {
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
                                   hierarchy->getPatchLevel(ln);
            tbox::Array<int> proc_mapping = 
               patch_level->getProcessorMapping().getProcessorMapping();

            for (pn = 0; pn < patch_level->getNumberOfPatches(); pn++) {
               int proc_num = proc_mapping[pn];
               int global_patch_id = getGlobalPatchNumber(hierarchy,ln,-1,pn);
               pms[global_patch_id].processor_number = proc_num;
               pms[global_patch_id].file_cluster_number = 
                  d_processor_in_file_cluster_number[proc_num];
               pms[global_patch_id].level_number = ln;
               pms[global_patch_id].patch_number = pn;
            }
         }
      }

      key_string = "patch_map";
      HDFputPatchMapStructArray(key_string,
                          pms,
                          tot_number_of_patches,
                          extents_group_id);

      delete [] pms;

      /*
       * Create "patch_extents" sub-database of extents group.
       */
      patchExtentsStruct *pes = new patchExtentsStruct[tot_number_of_patches];

      for (pn = 0; pn < tot_number_of_patches; pn++) {
         for (i = 0; i < VISIT_FIXED_DIM; i++) {
            pes[pn].lower[i] = 0;
            pes[pn].upper[i] = 0;
            pes[pn].xlo[i] = 0.;
            pes[pn].xhi[i] = 0.;
         }
      }


      /*
       * Set patch extents
       */
      if (d_grid_type != VISIT_DEFORMED) {
         //This is never entered in multiblock case
         tbox::Pointer< hier::PatchHierarchy<DIM> > patch_hierarchy =
            hierarchy;
         const tbox::Pointer< geom::CartesianGridGeometry<DIM> > ggeom = 
            patch_hierarchy->getGridGeometry();
         for (i = 0; i < DIM; i++) {
            geom_lo[i] = ggeom->getXLower()[i];
            dx_curr_lev[i] = ggeom->getDx()[i]; // coarsest level dx
         }

      } else {

         /*
          * Deformed grid - set extents to 0.
          */
         for (i = 0; i < DIM; i++) {
            geom_lo[i] = 0.;
            dx_curr_lev[i] = 0.; 
         }
      }

      if (d_is_multiblock) {

         int unique_level_number = 0;
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
               hierarchy->getPatchLevel(ln);
            for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
               tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
                  mblk_level->getPatchLevelForBlock(b);
               if (!patch_level.isNull()) {
                  const hier::BoxArray<DIM>& boxes = patch_level->getBoxes();

                  /*
                   * Set the dx for the next level
                   */
                  for (i = 0; i < DIM; i++) {
                     double scale_ratio =
                        (double)d_scaling_ratios[unique_level_number](i);
                     dx_curr_lev[i] = dx_curr_lev[i] / scale_ratio;
                  }

                  for (pn = 0; pn < boxes.getNumberOfBoxes(); pn++) {
                     int global_patch_id =
                        getGlobalPatchNumber(hierarchy,ln,b,pn);
                     const hier::Box<DIM>& box = boxes.getBox(pn);
                     const int* lower = box.lower();
                     const int* upper = box.upper();

                     for (i = 0; i < DIM; i++) {
                        pes[global_patch_id].lower[i] = lower[i];
                        pes[global_patch_id].upper[i] = upper[i];

                        patch_xlo = 0.;
                        patch_xhi = 0.;

                        pes[global_patch_id].xlo[i] = patch_xlo;
                        pes[global_patch_id].xhi[i] = patch_xhi;
                     }
                  }
                  unique_level_number++;
               }
            } // loop over patch boxes
         } // loop over levels
      } else {
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
               hierarchy->getPatchLevel(ln);
            const hier::BoxArray<DIM>& boxes = patch_level->getBoxes();

            /*
             * Set the dx for the next level
             */
            for (i = 0; i < DIM; i++) {
               double scale_ratio = (double)d_scaling_ratios[ln](i);
               dx_curr_lev[i] = dx_curr_lev[i] / scale_ratio;
            }

            for (pn = 0; pn < boxes.getNumberOfBoxes(); pn++) {
               int global_patch_id = getGlobalPatchNumber(hierarchy,ln,-1,pn);
               const hier::Box<DIM>& box = boxes.getBox(pn);
               const int* lower = box.lower();
               const int* upper = box.upper();

               for (i = 0; i < DIM; i++) {
                  pes[global_patch_id].lower[i] = lower[i];
                  pes[global_patch_id].upper[i] = upper[i];

                  if (d_grid_type != VISIT_DEFORMED) {
                     patch_xlo = geom_lo[i] + 
                        dx_curr_lev[i] * (double)lower[i];
                     patch_xhi = geom_lo[i] + 
                        dx_curr_lev[i] * (double)(upper[i]+1);
                  } else {
                     patch_xlo = 0.;
                     patch_xhi = 0.;
                  }
                  pes[global_patch_id].xlo[i] = patch_xlo;
                  pes[global_patch_id].xhi[i] = patch_xhi;   

               }
            } // loop over patch boxes
         } // loop over levels
      }

      /*
       * Write patch min/max for each variable.
       */
      tbox::Pointer<tbox::HDFDatabase> extents_materials_HDFGroup;
      for (typename tbox::List<VisItItem>::Iterator 
              ipi(d_plot_items); ipi; ipi++) {
         for (int comp = 0; comp < ipi().d_depth; comp++) {

            /*
             * Regular (i.e. not materials or species) variables
             */
            if (!(ipi().d_isa_material) && !(ipi().d_isa_species)) {

               key_string = ipi().d_visit_var_name[comp] + "-Extents";
               HDFputPatchMinMaxStructArray(
                  key_string,
                  ipi().d_master_min_max[comp],
                  tot_number_of_patches,
                  extents_group_id);

            } else if (ipi().d_isa_material) {

               /*
                * Create materials HDF group
                * set up HDF extents structure for materials:
                *  <extents group>
                *      materials group  (i.e. one constructed below)
                *         material_name group
                *            species group
                */
               if (extents_materials_HDFGroup.isNull()) {
                  extents_materials_HDFGroup = 
                     extents_HDFGroup->putDatabase("materials"); 
               }

               string mname = ipi().d_material_name;
               // material_name group
               tbox::Pointer<tbox::HDFDatabase> 
                  extents_material_name_HDFGroup =
                  extents_materials_HDFGroup->putDatabase(mname);
               hid_t extents_material_name_group_id = 
                      extents_material_name_HDFGroup->getGroupId();

               key_string = ipi().d_material_name + "-Fractions";
               HDFputPatchMinMaxStructArray(
                      key_string,
                      ipi().d_master_min_max[comp],
                      tot_number_of_patches,
                      extents_material_name_group_id);

               // species group
               if (ipi().d_species_names.getSize() > 0) {
                  ipi().d_extents_species_HDFGroup = 
                     extents_material_name_HDFGroup->putDatabase("species");
               }

            } else if (ipi().d_isa_species) {

               // species
               key_string = ipi().d_species_name;
               hid_t species_group_id = 
                  ipi().d_parent_material_pointer->
                  d_extents_species_HDFGroup->getGroupId(); // species group

               HDFputPatchMinMaxStructArray(
                  key_string,
                  ipi().d_master_min_max[comp],
                  tot_number_of_patches,
                  species_group_id);
            }

         } // loop over components
      } // loop over variables

      delete [] d_worker_min_max;
      
      key_string = "patch_extents";
      HDFputPatchExtentsStructArray(key_string,
                          pes,
                          tot_number_of_patches,
                          extents_group_id);
  
      delete [] pes;

      summary_HDFFilePointer->unmount();
 
   } // if VISIT_MASTER
   
   tbox::MPI::barrier();

   if (my_proc == VISIT_MASTER) {

      /*
       * Add this dump entry to dumps.visit file
       */
      if (d_time_step_number == 0) s_summary_file_opened = false;
      string path = d_top_level_directory_name + "/dumps.visit";
      string file = d_current_dump_directory_name + "/" + d_summary_filename;

      /*
       * If summary file has not yet been opened, open and write file.
       * If it has been opened, append this to file.
       */
      if (!s_summary_file_opened) {
         s_summary_file_opened = true;
         ofstream sfile(path.c_str(), ios::out);
         sfile <<  file << "\n";
         sfile.close(); 
      } else {
         ofstream sfile(path.c_str(), ios::app);
         sfile <<  file << "\n";
         sfile.close(); 
      }
   }

   tbox::MPI::barrier();
}

/*
*************************************************************************
*                                                                       *
* Private function to store min/max information on each patch for each  *
* variable.  The "master" processor allocates an array that will hold   *
* the global data.  The "worker" processors send this information to    *
* the master which in turn unpacks and stores the data.                 *
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::exchangeMinMaxPatchInformation(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int coarsest_plot_level,
   const int finest_plot_level)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(coarsest_plot_level>=0);
   assert(finest_plot_level>=0);
#endif

   /*
    * Compute max number of patches on any processor, and the total number of
    * patches in the problem.
    */
   int ln, pn, comp, item_ctr;
   int number_local_patches = 0;
   int tot_number_of_patches = 0;

   if (d_is_multiblock) {
      for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
            hierarchy->getPatchLevel(ln);
         for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               tot_number_of_patches += patch_level->getNumberOfPatches();
               for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
                    ip; ip++) {
                  number_local_patches++;
               }
            }
         }
      }
   } else {
      for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
            hierarchy->getPatchLevel(ln);
         tot_number_of_patches += patch_level->getNumberOfPatches();
         for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
              ip; ip++) {
            number_local_patches++;
         }
      }
   }

   int max_number_local_patches = tbox::MPI::maxReduction(number_local_patches);
   
   int num_items_to_plot = d_number_visit_variables_plus_depth
      + d_materials_names.getSize() // number materials 
      + d_number_species;
   
   int message_size = max_number_local_patches * 
      sizeof(patchMinMaxStruct) * num_items_to_plot;

   if (tbox::MPI::getRank() != VISIT_MASTER) {
      
      /*
       * Worker processor:  send contents of d_worker_min_max array that
       * was setup in "initializePlotVariableMinMaxInfo()" to the
       * master processor.
       */
      if (number_local_patches > 0) {
         tbox::MPI::sendBytes((void *)d_worker_min_max, 
                             message_size,
                             VISIT_MASTER);
      }
      
   } else { // (tbox::MPI::getRank() == VISIT_MASTER)
      
      /*
       * Master processor:  Receive the min/max information sent by the 
       * worker processors (above).  Unpack into the local d_mm array
       * for each plot variable.
       */

      // recv buffer large enough to receive info from any processor. 
      patchMinMaxStruct *buf = NULL;
      if (d_number_working_slaves > 0) {
         buf = new patchMinMaxStruct[max_number_local_patches 
                                     * num_items_to_plot]; 
      }
      
      /*
       * Receive information sent by "sending_proc".
       */
      int number_msgs_recvd = 0;
      while (number_msgs_recvd < d_number_working_slaves) {
         int sending_proc = tbox::MPI::recvBytes((void *)buf,
                                                message_size);
         number_msgs_recvd++;
         
         /*
          * Unpack the information from buf and fill d_worker_min_max array.
          */ 
         item_ctr = 0;

         if (d_is_multiblock) {
            for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
               tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
                  hierarchy->getPatchLevel(ln);
               for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
                  tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
                     mblk_level->getPatchLevelForBlock(b);
                  if (!patch_level.isNull()) {
                     tbox::Array<int> proc_mapping =
                        patch_level->
                           getProcessorMapping().getProcessorMapping();

                     int npatches_on_level = proc_mapping.getSize();
                     for (pn = 0; pn < npatches_on_level; pn++) {
                        if (proc_mapping[pn] == sending_proc) {
                           int global_patch_id =
                              getGlobalPatchNumber(hierarchy, ln, b, pn);
                           for (typename tbox::List<VisItItem>::Iterator
                                   ipi(d_plot_items); ipi; ipi++) {
                              for (comp = 0; comp < ipi().d_depth; comp++) {
                                 ipi().d_master_min_max[comp][global_patch_id] =
                                    buf[item_ctr];
                                 item_ctr++;
                              }
                           }  // variables
                        } // patch from sending proc?
                     } // patches
                  }
               }
            } // levels
         } else {
            for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
               tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
                  hierarchy->getPatchLevel(ln);
               tbox::Array<int> proc_mapping = 
                  patch_level->getProcessorMapping().getProcessorMapping();

               int npatches_on_level = proc_mapping.getSize();
               for (pn = 0; pn < npatches_on_level; pn++) {
                  if (proc_mapping[pn] == sending_proc) {
                     int global_patch_id = 
                        getGlobalPatchNumber(hierarchy, ln,-1,  pn);
                     for (typename tbox::List<VisItItem>::Iterator 
                             ipi(d_plot_items); ipi; ipi++) {
                        for (comp = 0; comp < ipi().d_depth; comp++) {
                           ipi().d_master_min_max[comp][global_patch_id] = 
                              buf[item_ctr];
                           item_ctr++;
                        }
                     }  // variables
                  } // patch from sending proc?
               } // patches
            } // levels
         }
      } // while msgs_recvd < working procs

      if (d_number_working_slaves > 0) {
         delete [] buf;
      }

   } // proc == VISIT_MASTER?
   
}   


/*
*************************************************************************
*                                                                       *
* Private function to find and write to summary file parent & child     *
* info.  For each global patch number, find its children using the      *
* box_tree.  Record the children, as well as each child's parent, in a  *
* child_parent array. A child_ptrs array records for each global patch  *
* number, the number of children that patch has, as well as the offset  *
* into the child_parent array where the patch numbers of those children *
* are stored.  If a patch has no children, offset = -1. Next, the child *
* info from the child_parent array is copied into the final child array.*  
* Then the child_parent array is sorted by child number.  Now all       *
* parents of a given patch are grouped together in this sorted array.   *
* The parents are then stored in a parent array, and a parent_ptrs      *
* array is created similar to the child_ptrs array.                     *
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::writeParentChildInfoToSummaryHDFFile(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   tbox::Pointer<tbox::HDFDatabase> basic_HDFGroup)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
#endif

   struct cpPointerStruct {  // auxiliary info for child/parent data struct
      int offset;
      union {
         int number_children;
         int number_parents;
      }u;
   };

   /*
    * Find child patches of each global patch number
    */
   int tot_number_of_patches = 0;
   int finest_level = hierarchy->getFinestLevelNumber();
   int ln;

   if (d_is_multiblock) {
      for (ln = 0; ln <= finest_level; ln++) {
         tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
            hierarchy->getPatchLevel(ln);
         for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mblk_level->getPatchLevelForBlock(b);
            if (!patch_level.isNull()) {
               tot_number_of_patches += patch_level->getNumberOfPatches();
            }
         }
      }
   } else {
      for (ln = 0; ln <= finest_level; ln++) {
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
            hierarchy->getPatchLevel(ln);
         tot_number_of_patches += patch_level->getNumberOfPatches();
      }
   }
   int chunk_size = 2 * tot_number_of_patches;
   int current_child_parent_max_size = chunk_size;
   struct cpPointerStruct *child_ptrs = 
      new struct cpPointerStruct[tot_number_of_patches];
   struct childParentStruct *child_parent = 
      new struct childParentStruct[chunk_size];
   int child_parent_idx = 0;
   int child_ptrs_idx = 0;

   if (d_is_multiblock) {
      tbox::Pointer< mblk::MultiblockPatchLevel<DIM> > mblk_level =
         hierarchy->getPatchLevel(ln);
      for (int b = 0; b < mblk_level->getNumberBlocks(); b++) {
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
            mblk_level->getPatchLevelForBlock(b);
         if (!patch_level.isNull()) {
            hier::BoxArray<DIM> coarser_boxes = patch_level->getBoxes();

            const int num_coarser_boxes = coarser_boxes.getNumberOfBoxes();
            tbox::Pointer< hier::BoxTree<DIM> > child_box_tree;
            bool finest_level_in_block = false;

            if (ln != finest_level) {
               tbox::Pointer< mblk::MultiblockPatchLevel<DIM> >
                  child_mblk_level =
                     hierarchy->getPatchLevel(ln + 1);
               tbox::Pointer< hier::PatchLevel<DIM> > child_patch_level =
                  child_mblk_level->getPatchLevelForBlock(b);

               if (!child_patch_level.isNull()) {
                  coarser_boxes.refine(
                     child_patch_level->getRatioToCoarserLevel());
                  child_box_tree = child_patch_level->getBoxTree();
               } else {
                  finest_level_in_block = true;
               }
            }

            for (int icp = 0; icp < num_coarser_boxes; icp++) {
               if (ln == finest_level || finest_level_in_block) {
                  child_ptrs[child_ptrs_idx].u.number_children = 0;
                  child_ptrs[child_ptrs_idx++].offset =
                     VISIT_UNDEFINED_INDEX;
               } else {
                  tbox::Array<int> child_patch_array;
                  child_box_tree->findOverlapIndices(
                     child_patch_array,
                     coarser_boxes.getBox(icp));
                  int num_kids = child_patch_array.getSize();
                  child_ptrs[child_ptrs_idx].u.number_children = num_kids;

                  if (num_kids == 0) {
                     child_ptrs[child_ptrs_idx++].offset =
                        VISIT_UNDEFINED_INDEX;
                  }
                  else {
                     child_ptrs[child_ptrs_idx++].offset = child_parent_idx;
                     if ((child_parent_idx + num_kids) >
                                        current_child_parent_max_size) {
                        current_child_parent_max_size += chunk_size;
                        struct childParentStruct *temp = child_parent;
                        child_parent =
                            new struct
                               childParentStruct[current_child_parent_max_size];
                        for (int idx = 0; idx < child_parent_idx; idx++) {
                           child_parent[idx].child = temp[idx].child;
                           child_parent[idx].parent = temp[idx].parent;
                        }
                        delete [] temp;
                     }

                     for (int idx = 0; idx < num_kids; idx++) {
                        child_parent[child_parent_idx].child =
                           getGlobalPatchNumber(hierarchy, ln + 1,
                                                b, child_patch_array[idx]);
                        child_parent[child_parent_idx++].parent =
                           getGlobalPatchNumber(hierarchy, ln, b, icp);
                     }
                  }
               }
            }
         }
      }
   } else {
      for (ln = 0; ln <= finest_level; ln++) {
         tbox::Pointer<hier::PatchLevel<DIM> > patch_level =
            hierarchy->getPatchLevel(ln);
         hier::BoxArray<DIM> coarser_boxes = patch_level->getBoxes();

         const int num_coarser_boxes = coarser_boxes.getNumberOfBoxes();
         tbox::Pointer< hier::BoxTree<DIM> > child_box_tree;

         if (ln != finest_level) {
            tbox::Pointer< hier::PatchLevel<DIM> > child_patch_level = 
               hierarchy->getPatchLevel(ln + 1);
            coarser_boxes.refine(child_patch_level->getRatioToCoarserLevel());
            child_box_tree = child_patch_level->getBoxTree();
         }
      
         for (int icp = 0; icp < num_coarser_boxes; icp++) {
            if (ln == finest_level) {
               child_ptrs[child_ptrs_idx].u.number_children = 0;
               child_ptrs[child_ptrs_idx++].offset = VISIT_UNDEFINED_INDEX;
            } else {
               tbox::Array<int> child_patch_array;
               child_box_tree->findOverlapIndices(child_patch_array, 
                                                  coarser_boxes.getBox(icp));
               int num_kids = child_patch_array.getSize();
               child_ptrs[child_ptrs_idx].u.number_children = num_kids;

               if (num_kids == 0) {
                  child_ptrs[child_ptrs_idx++].offset = VISIT_UNDEFINED_INDEX;
               }
               else {
                  child_ptrs[child_ptrs_idx++].offset = child_parent_idx;
                  if ((child_parent_idx + num_kids) > 
                                     current_child_parent_max_size) {
                     current_child_parent_max_size += chunk_size;
                     struct childParentStruct *temp = child_parent;
                     child_parent = 
                         new struct 
                            childParentStruct[current_child_parent_max_size];
                     for (int idx = 0; idx < child_parent_idx; idx++) {
                        child_parent[idx].child = temp[idx].child;
                        child_parent[idx].parent = temp[idx].parent;
                     }
                     delete [] temp;
                  }

                  for (int idx = 0; idx < num_kids; idx++) {
                     child_parent[child_parent_idx].child = 
                        getGlobalPatchNumber(hierarchy, ln + 1,
                                             -1, child_patch_array[idx]);
                     child_parent[child_parent_idx++].parent = 
                        getGlobalPatchNumber(hierarchy, ln, -1, icp);
                  }
               }
            }
         }
      }
   }

   int *parent_array = (int*)NULL;
   int *child_array = (int*)NULL; 
   int parent_array_length = 0;
   int child_array_length = child_parent_idx;

   struct cpPointerStruct *parent_ptrs = (cpPointerStruct*)NULL;

   // copy child info to child array
   if (child_array_length > 0) {
      child_array = new int[child_array_length];
      for (int idx = 0; idx < child_array_length; idx++) {
         child_array[idx] = child_parent[idx].child;
      }

      // sort child_parent array by child patch number.
      qsort((char *)child_parent,
            child_array_length,
            sizeof(struct childParentStruct), 
            &childParentCompareFunc);

      // now record parents in the parent array
      parent_array = new int[chunk_size];
      int parent_size = chunk_size;
      int cp_idx = 0;
      int next_parent = 0;
      parent_ptrs = 
         new struct cpPointerStruct[tot_number_of_patches];

      for (int gpn = 0; gpn < tot_number_of_patches; gpn++) {
         if (gpn < child_parent[cp_idx].child) {
            parent_ptrs[gpn].offset = VISIT_UNDEFINED_INDEX;
            parent_ptrs[gpn].u.number_parents = 0;
         } else {
            int num_pars = 0;
            while (child_parent[cp_idx].child == gpn) {
               if (num_pars == 0) {
                  parent_ptrs[gpn].offset = next_parent;
               }
               num_pars++;
               if (next_parent >= parent_size) {
                  // increase size of parent_array
                  int old_parent_size = parent_size;
                  int *temp = new int[old_parent_size];
                  for (int i = 0; i < old_parent_size; i++) {
                     temp[i] = parent_array[i];
                  }
                  delete [] parent_array;
                  parent_size += chunk_size;
                  parent_array = new int [parent_size];
                  for (int i = 0; i < old_parent_size; i++) {
                     parent_array[i] = temp[i];
                  }
                  delete [] temp;
               }
               parent_array[next_parent++] = child_parent[cp_idx++].parent;
            }
            parent_ptrs[gpn].u.number_parents = num_pars;
         }
      }
      parent_array_length = next_parent;
   }
   delete [] child_parent;

   /*
    * Write parent & child info to summary file
    */
   string key_string = "child_array_length";
   basic_HDFGroup->putInteger(key_string, child_array_length);
   key_string = "parent_array_length";
   basic_HDFGroup->putInteger(key_string, parent_array_length);
   hid_t basic_group_id = basic_HDFGroup->getGroupId();
   if (child_array_length > 0) {
      key_string = "child_array";
      basic_HDFGroup->putIntegerArray(key_string,
                                      child_array,
                                      child_array_length);
      
      key_string = "child_pointer_array";
      HDFputChildParentStructArray(
                   key_string,
                   (void *) child_ptrs,
                   tot_number_of_patches,
                   basic_group_id,
                   sizeof(cpPointerStruct),
                   "number_children");
   }
   if (child_array) delete [] child_array;
   if (child_ptrs) delete [] child_ptrs;

   if (parent_array_length > 0) {
      key_string = "parent_array";
      basic_HDFGroup->putIntegerArray(key_string,
                                      parent_array,
                                      parent_array_length);
      key_string = "parent_pointer_array";
      HDFputChildParentStructArray(
                   key_string,
                   (void *) parent_ptrs,
                   tot_number_of_patches,
                   basic_group_id,
                   sizeof(cpPointerStruct),
                   "number_parents");
   }

   if (parent_array) delete [] parent_array;
   if (parent_ptrs) delete [] parent_ptrs;
   
}


/*
*************************************************************************
*                                                                       *
*    childParentCompareFunc() used by qsort to sort child_parent array  *
*        by child patch num to find all parents of a given child        *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VisItDataWriter<DIM>::childParentCompareFunc(
   const void *s1, 
   const void *s2)
{
   struct childParentStruct* a = (struct childParentStruct*) s1;
   struct childParentStruct* b = (struct childParentStruct*) s2;
   if (((struct childParentStruct*)a)->child > 
       ((struct childParentStruct*)b)->child) {
      return (1);
   } else if (((struct childParentStruct*)b)->child > 
              ((struct childParentStruct*)a)->child) {
      return (VISIT_UNDEFINED_INDEX);
   } else {
      return(0);
   }
}

/*
*************************************************************************
*                                                                       *
* Private utility function to pack DIM patch data into 1D double       *
* precision buffer, omitting ghost data if necessary. Data is packed in *
* column major order, i.e. x0,y0,z0, x1,y0,z0, x2,y0,z0, ...            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::packPatchDataIntoDoubleBuffer(
   const tbox::Pointer< hier::PatchData<DIM> > pdata,
   const int depth_index,
   const variable_data_type data_type,
   const hier::Box<DIM> patch_box,
   double* buffer,
   const variable_centering centering)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth_index>=0);
   assert((patch_box * pdata->getGhostBox()) == patch_box);
#endif

   /*
    * Currently, Fortran calls in this method limit it to two or 
    * three dimensions only.  It will not work for DIM < 2 or 
    * DIM > 3.  Give the user a run-time error asserting this.
    */
   if (DIM < 2 || DIM > 3) {
      TBOX_ERROR(d_object_name << ":packPatchDataIntoDoubleBuffer()"
		 << "\n  This case has DIM = " << DIM
		 << "\n  Dimensions < 2 or > 3 are not supported at"
		 << "\n  this time." << endl);   
   }

   int buf_size = getBufferSize(patch_box,
                                hier::IntVector<DIM>(0),
                                centering);

   hier::Index<DIM> databox_lower = pdata->getGhostBox().lower();
   hier::Index<DIM> databox_upper = pdata->getGhostBox().upper();

   hier::Box<DIM> plot_box = patch_box;
   hier::Index<DIM> plolower = plot_box.lower();
   hier::Index<DIM> ploupper = plot_box.upper();

   /*
    * Extend index boundaries by one point if NODE data is used.
    */
   if (centering == VISIT_NODE || centering == VISIT_UNKNOWN_NODE) {
      databox_upper += 1;
      ploupper += 1;
   }
   
   switch(data_type) {

#ifdef HAVE_FLOAT
      case VISIT_FLOAT: {
         float *dat_ptr = NULL;
         if (centering == VISIT_CELL) {
            const tbox::Pointer< pdat::CellData<DIM,float> > fpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!fpdata.isNull());
#endif
            dat_ptr = fpdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<DIM,float> cell_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         }
         else if (centering == VISIT_NODE) {
            const tbox::Pointer< pdat::NodeData<DIM,float> > fpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!fpdata.isNull());
#endif
            dat_ptr = fpdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<DIM,float> node_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
	 if (DIM == 2) {
	    cpfdat2buf2d_(databox_lower(0), databox_lower(1),
			  plolower(0), plolower(1),
			  ploupper(0), ploupper(1),
			  databox_upper(0), databox_upper(1),
			  dat_ptr, buffer, buf_size);
	 } else if (DIM == 3) {
	    cpfdat2buf3d_(databox_lower(0), databox_lower(1), databox_lower(2),
			  plolower(0), plolower(1), plolower(2),
			  ploupper(0), ploupper(1), ploupper(2),
			  databox_upper(0), databox_upper(1), databox_upper(2),
			  dat_ptr, buffer, buf_size);
	 }
         break;
      } 
#endif  
         
      case VISIT_DOUBLE: {
         double *dat_ptr = NULL;
         if (centering == VISIT_CELL) {
            const tbox::Pointer< pdat::CellData<DIM,double> > dpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!dpdata.isNull());
#endif
            dat_ptr = dpdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<DIM,double> cell_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         }
         else if (centering == VISIT_NODE) {
            const tbox::Pointer< pdat::NodeData<DIM,double> > dpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!dpdata.isNull());
#endif
            dat_ptr = dpdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<DIM,double> node_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
	 if (DIM == 2) {
	    cpddat2buf2d_(databox_lower(0), databox_lower(1),
			  plolower(0), plolower(1),
			  ploupper(0), ploupper(1),
			  databox_upper(0), databox_upper(1),
			  dat_ptr, buffer, buf_size); 
	 } else if (DIM == 3) {
	    cpddat2buf3d_(databox_lower(0), databox_lower(1), databox_lower(2),
			  plolower(0), plolower(1), plolower(2),
			  ploupper(0), ploupper(1), ploupper(2),
			  databox_upper(0), databox_upper(1), databox_upper(2),
			  dat_ptr, buffer, buf_size);
	 }
         break;
      }

      case VISIT_INT: {
         int *dat_ptr = NULL;
         if (centering == VISIT_CELL) {
            const tbox::Pointer< pdat::CellData<DIM,int> > ipdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!ipdata.isNull());
#endif
            dat_ptr = ipdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<DIM,int> cell_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         }
         else if (centering == VISIT_NODE) {
            const tbox::Pointer< pdat::NodeData<DIM,int> > ipdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!ipdata.isNull());
#endif
            dat_ptr = ipdata->getPointer(depth_index);
         }
         else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<DIM,int> node_copy(plot_box, 1, hier::IntVector<DIM>(0));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
	 if (DIM == 2) {
	    cpidat2buf2d_(databox_lower(0), databox_lower(1),
			  plolower(0), plolower(1),
			  ploupper(0), ploupper(1),
			  databox_upper(0), databox_upper(1),
			  dat_ptr, buffer, buf_size); 
	 } else if (DIM == 3) {
	    cpidat2buf3d_(databox_lower(0), databox_lower(1), databox_lower(2),
			  plolower(0), plolower(1), plolower(2),
			  ploupper(0), ploupper(1), ploupper(2),
			  databox_upper(0), databox_upper(1), databox_upper(2),
			  dat_ptr, buffer, buf_size);
	 }
         break;
      }

      default: {
         TBOX_ERROR(d_object_name << ":packPatchDataIntoDoubleBuffer()"
                    << "\n  Unknown type.  ***Exiting." << endl);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Create a 2D integer array entry in an HDF database with the specified *
* key name.  The array type is based on the hdf type H5T_NATIVE_INT.    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::HDFputIntegerArray2D(
   const string& key, 
   const int *data, 
   const int nelements0,
   const int nelements1,
   const hid_t group_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert((nelements0 > 0) && (nelements1 > 0));
#endif
   herr_t errf;
   if ((nelements0 > 0) && (nelements1 > 0)) {
      hsize_t dim[] = {nelements0, nelements1};
      hid_t space = H5Screate_simple(2, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                H5T_NATIVE_INT, 
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, 
                      H5T_NATIVE_INT, 
                      H5S_ALL, 
                      H5S_ALL, 
                      H5P_DEFAULT, 
                      data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputIntegerArray2D()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a 2D double array entry in an HDF database with the specified  *
* key name.  The array type is based on the hdf type H5T_NATIVE_DOUBLE. *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::HDFputDoubleArray2D(
   const string& key, 
   const double *data, 
   const int nelements0,
   const int nelements1,
   const hid_t group_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert((nelements0 > 0) && (nelements1 > 0));
#endif
   herr_t errf;
   if ((nelements0 > 0) && (nelements1 > 0)) {
      hsize_t dim[] = {nelements0, nelements1};
      hid_t space = H5Screate_simple(2, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                H5T_NATIVE_DOUBLE,
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, 
                      H5T_NATIVE_DOUBLE, 
                      H5S_ALL, 
                      H5S_ALL, 
                      H5P_DEFAULT, 
                      data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputDoubleArray2D()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}



/*
*************************************************************************
*                                                                       *
* Create an array of patch extent (pe) structs in an HDF database       *
* with the specified key name.                                          *
*                                                                       *
*************************************************************************
*/


template<int DIM> void VisItDataWriter<DIM>::HDFputPatchExtentsStructArray(
   const string& key, 
   const patchExtentsStruct *data,
   const int nelements,
   const hid_t group_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert(nelements > 0);
#endif
   herr_t errf;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t pe_id = H5Tcreate(H5T_COMPOUND, sizeof(patchExtentsStruct));
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(pe_id >= 0);
#endif
      hsize_t dim1[1];
      dim1[0] = VISIT_FIXED_DIM;  
      hid_t intXdType = H5Tarray_create(H5T_NATIVE_INT, 1, dim1, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(intXdType >= 0);
#endif
      hid_t doubleXdType = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim1, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(doubleXdType >= 0);
#endif
      errf = H5Tinsert(pe_id, 
                       "lower",
                       HOFFSET(patchExtentsStruct, lower),
                       intXdType);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pe_id, 
                       "upper",
                       HOFFSET(patchExtentsStruct, upper),
                       intXdType);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pe_id, 
                       "xlo",
                       HOFFSET(patchExtentsStruct, xlo),
                       doubleXdType); 
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pe_id, 
                       "xup",
                       HOFFSET(patchExtentsStruct, xhi),
                       doubleXdType); 
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                pe_id, 
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, 
                      pe_id, 
                      H5S_ALL, 
                      H5S_ALL,
                      H5P_DEFAULT, 
                      data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(pe_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(intXdType);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(doubleXdType);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif


   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputPatchExtentsStructArray()"
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create patch map for each of the patch extents HDF entries            *
* with the specified key name.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::HDFputPatchMapStructArray(
   const string& key, 
   const patchMapStruct *data,
   const int nelements,
   const hid_t group_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert(nelements > 0);
#endif
   herr_t errf;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t pm_id = H5Tcreate(H5T_COMPOUND, sizeof(patchMapStruct));
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(pm_id >= 0);
#endif
      errf = H5Tinsert(pm_id, 
                       "processor_number",
                       HOFFSET(patchMapStruct, processor_number),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pm_id, 
                       "file_cluster_number",
                       HOFFSET(patchMapStruct, file_cluster_number),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pm_id, 
                       "level_number",
                       HOFFSET(patchMapStruct, level_number),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS      
      assert(errf >= 0);
#endif
      errf = H5Tinsert(pm_id, 
                       "patch_number",
                       HOFFSET(patchMapStruct, patch_number),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                pm_id, 
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, 
                      pm_id, 
                      H5S_ALL, 
                      H5S_ALL,
                      H5P_DEFAULT, 
                      data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(pm_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif

   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputPatchMapStructArray()"
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create an array of max-min double (mm) structs an HDF database with   *
* the specified key name.                                               *
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::HDFputPatchMinMaxStructArray(
   const string& key, 
   const patchMinMaxStruct *data,
   const int nelements,
   const hid_t group_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert(nelements > 0);
#endif
   herr_t errf;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(patchMinMaxStruct));
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(s1_tid >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "data_is_defined",
                       HOFFSET(patchMinMaxStruct, patch_data_on_disk),
                       H5T_NATIVE_CHAR);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "material_composition_flag",
                       HOFFSET(patchMinMaxStruct, material_composition_code),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "species_composition_flag",
                       HOFFSET(patchMinMaxStruct, species_composition_code),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "min",
                       HOFFSET(patchMinMaxStruct, min),
                       H5T_NATIVE_DOUBLE);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "max",
                       HOFFSET(patchMinMaxStruct, max),
                       H5T_NATIVE_DOUBLE);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                s1_tid, 
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, 
                      s1_tid, 
                      H5S_ALL, 
                      H5S_ALL,
                      H5P_DEFAULT, 
                      data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(s1_tid);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif

   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputPatchMinMaxStructArray()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}


/*
*************************************************************************
*                                                                       *
* Create an array of child/parent pointer (cpp) structs an HDF          *
* database with the specified key name.                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VisItDataWriter<DIM>::HDFputChildParentStructArray(
   const string& key, 
   const void *data, 
   const int nelements,
   hid_t group_id,
   const int sizeOfStruct,
   const string& field_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != NULL);
   assert(nelements > 0);
#endif
   herr_t errf;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(space >= 0);
#endif
      hid_t s1_tid = H5Tcreate(H5T_COMPOUND, sizeOfStruct);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(s1_tid >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       "offset",
                       0,
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tinsert(s1_tid, 
                       field_name.c_str(),
                       sizeof(int),
                       H5T_NATIVE_INT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      hid_t dataset = H5Dcreate(group_id, 
                                key.c_str(), 
                                s1_tid, 
                                space, 
                                H5P_DEFAULT);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Sclose(space);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Tclose(s1_tid);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif
      errf = H5Dclose(dataset);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(errf >= 0);
#endif


   } else {
      TBOX_ERROR("VisItDataWriter<DIM>::HDFputChildParentStructArray()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Private function to calculate required size for a buffer              *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VisItDataWriter<DIM>::getBufferSize(
   const hier::Box<DIM> patch_box, 
   const hier::IntVector<DIM>& ghost_cell_width, 
   const variable_centering centering)
{
   int buf_size = 1;
   int cen = 0;
   if (centering == VISIT_CELL || centering == VISIT_UNKNOWN_CELL) {
      cen = 1;
   }
   else if (centering == VISIT_NODE || centering == VISIT_UNKNOWN_NODE) {
      cen = 2;
   }
   const int* lower = patch_box.lower();
   const int* upper = patch_box.upper();

   for (int i = 0; i < DIM; i++) {
      buf_size *= upper[i]-lower[i]+cen+(2*ghost_cell_width(i));
   }
   return buf_size;
}


/*
*************************************************************************
*                                                                       *
* Dump plotitem fields for debugging purposes                           *
*                                                                       *
*************************************************************************
*/
template<int DIM> void VisItDataWriter<DIM>::dumpItem(
   VisItItem& plotitem,
   ostream& os) const
{
   os << "d_var_name: " << plotitem.d_var_name << "\n";
   string type;
   if (plotitem.d_var_type == VISIT_SCALAR) {
      type = "SCALAR";
   } else if (plotitem.d_var_type == VISIT_VECTOR) {
      type = "VECTOR";
   } else if (plotitem.d_var_type == VISIT_TENSOR) {
      type = "TENSOR";
   }
   os << "d_var_type: " << type << "\n";
   string data_type;
   if (plotitem.d_var_data_type == VISIT_DOUBLE) {
      data_type = "DOUBLE";
   } else if (plotitem.d_var_data_type == VISIT_FLOAT) {
      data_type = "FLOAT";
   } else if (plotitem.d_var_data_type == VISIT_INT) {
      data_type = "INT";
   }
   os << "d_var_data_type: " << data_type << "\n";
   string cent;
   if (plotitem.d_var_centering == VISIT_CELL) {
      cent = "CELL";
   } else if (plotitem.d_var_centering == VISIT_UNKNOWN_CELL) {
      cent = "UNKNOWN_CELL";
   } else if (plotitem.d_var_centering == VISIT_NODE) {
      cent = "NODE";
   } else if (plotitem.d_var_centering == VISIT_UNKNOWN_NODE) {
      cent = "VISIT_UNKNOWN_NODE";
   }
   os << "d_var_centering: " << cent << "\n";
   os << "d_patch_data_index: " << plotitem.d_patch_data_index << "\n";
   os << "d_depth: " << plotitem.d_depth << "\n";
   os << "d_start_depth_index: " << plotitem.d_start_depth_index << "\n";
   os << "d_scale_factor: " << plotitem.d_scale_factor << "\n";
   os << "d_derived_writer ptr: " << plotitem.d_derived_writer << "\n";
   
   int i;
   for (i = 0; i < plotitem.d_depth; i++) {
      os << "   comp_name[" << i << "]: " 
           << plotitem.d_visit_var_name[i] << "\n";
   }

   os << "d_isa_material: " << plotitem.d_isa_material << "\n";
   os << "d_material_name: " << plotitem.d_material_name << "\n";
   os << "d_materials_writer ptr: " << plotitem.d_materials_writer << "\n";
   if (plotitem.d_isa_material) {
      os << "  Species for this material: " << "\n";
      for (i = 0; i < plotitem.d_species_names.getSize(); i++) {
         os << "   species[" << i << "]: " 
              << plotitem.d_species_names[i] << "\n";
      }
   }
   
   os << "d_isa_species: " << plotitem.d_isa_species << "\n";
   if (plotitem.d_parent_material_pointer) {
     os << "parent_material: " 
          << plotitem.d_parent_material_pointer->d_material_name << "\n";
   } else {
     os << "parent_material: Not Applicable\n";
   }
}

}
}
#endif

#endif
