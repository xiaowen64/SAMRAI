//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/apputils/plotting/CartesianVizamraiDataWriter.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1889 $
// Modified:    $LastChangedDate: 2008-01-22 16:46:52 -0800 (Tue, 22 Jan 2008) $
// Description: Simple tool to facilitate dumping data to file for Vizamrai
//

#ifndef included_appu_CartesianVizamraiDataWriter_C
#define included_appu_CartesianVizamraiDataWriter_C

#include "CartesianVizamraiDataWriter.h"

#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "PatchDescriptor.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "tbox/ArenaManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define VIZAMRAI_FILE_FORMAT_VERSION (6)

#define DUMP_DOUBLE_DATA   (0)
#define DUMP_FLOAT_DATA    (1)

#define DATA_TO_VIZ_IS_UNDEFINED (-1)
#define DATA_TO_VIZ_IS_DOUBLE    (0)
#define DATA_TO_VIZ_IS_FLOAT     (1)
#define DATA_TO_VIZ_IS_INT       (2)

#define DOMAIN_SCALE_FACTOR      (0.01)

namespace SAMRAI {
    namespace appu {

/*
*************************************************************************
*                                                                       *
* The constructor simply sets default object state.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CartesianVizamraiDataWriter<DIM>::CartesianVizamraiDataWriter(
   const std::string& name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif
   d_object_name = name;
   d_wrote_plot_file = false; 
   d_finest_plot_level_now = -1;
   d_finest_plot_level_in_all = -1;
   d_scaling_ratios.resizeArray(1);
   d_scaling_ratios[0] = hier::IntVector<DIM>(1);
   d_plot_type = DUMP_FLOAT_DATA;
   d_default_derived_writer = NULL;
   d_domain_scale_length = 0.0;
}

/*
*************************************************************************
*                                                                       *
* The destructor implicitly deallocates the list of plot data.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CartesianVizamraiDataWriter<DIM>::~CartesianVizamraiDataWriter()
{
}

/*
*************************************************************************
*                                                                       *
* Set maximum number of hierarchy levels to write to any plot file.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setFinestLevelToPlot(
   int finest_level_number)
{
   if ( d_wrote_plot_file && 
        (finest_level_number > d_scaling_ratios.getSize()+1) ) {
      TBOX_WARNING("CartesianVizamraiDataWriter<DIM>::setFinestLevelToPlot"
         << " warning...\n"
         << "    Data writer with name " << d_object_name
         << "\n     Attempting to change finest plot level from" 
         << d_finest_plot_level_in_all 
         << " to " << finest_level_number 
         << "\n     but scaling ratio data has not yet been given" 
         << "\n     Thus, finest level in any plot file will be level " 
         << d_finest_plot_level_in_all << std::endl);
   } else {
      if (!d_wrote_plot_file) {
         d_finest_plot_level_in_all = finest_level_number;
      }
      d_finest_plot_level_now = finest_level_number;
   }
}   

/*
*************************************************************************
*                                                                       *
* Set scaling ratio information for given level.                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setRatioToCoarserLevel(
   int level_number,
   const hier::IntVector<DIM>& ratio)
{
   if (d_wrote_plot_file) {
      TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::setRatioToCoarserLevel"
         << " error...\n"
         << "    Data writer with name " << d_object_name
         << "\n     Attempting to set level scaling information after"
         << "\n     writing a plot file.  This is not allowed." << std::endl);
   }

   if (level_number > d_scaling_ratios.getSize()-1) {
      d_scaling_ratios.resizeArray(level_number+1);
   }
   d_scaling_ratios[level_number] = ratio;
   d_finest_plot_level_in_all = 
      tbox::MathUtilities<int>::Max( d_finest_plot_level_in_all,
                                     level_number );
   d_finest_plot_level_now = d_finest_plot_level_in_all; 
}

/*
*************************************************************************
*                                                                       *
* Set format of plot data.                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setPlotDataToDouble()
{
   if (d_wrote_plot_file && (d_plot_type != DUMP_DOUBLE_DATA)) {
      TBOX_WARNING("CartesianVizamraiDataWriter<DIM>::setPlotDataToDouble"
         << " warning...\n"
         << "    data writer with name " << d_object_name
         << "\n     Changing plot type from float to double"
         << "\n     after writing plot data to file" << std::endl);
   }
   d_plot_type = DUMP_DOUBLE_DATA; 
}

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setPlotDataToFloat()
{
   if (d_wrote_plot_file && (d_plot_type != DUMP_FLOAT_DATA)) {
      TBOX_WARNING("CartesianVizamraiDataWriter<DIM>::setPlotDataToFloat"
         << " warning...\n"
         << "    data writer with name " << d_object_name
         << "\n     Changing plot type from double to float"
         << "\n     after writing plot data to file" << std::endl);
   }
   d_plot_type = DUMP_FLOAT_DATA; 
}

/*
*************************************************************************
*                                                                       *
* Set default user-defined derived data writer.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setDerivedDataWriter(
   VisDerivedDataStrategy<DIM>* derived_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(derived_writer != (VisDerivedDataStrategy<DIM>*)NULL);
#endif
   d_default_derived_writer = derived_writer;
}

/*
*************************************************************************
*                                                                       *
* Set default directory name for writing Vizamrai files.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::setDirectoryName(
   const std::string& directory_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!directory_name.empty());
#endif
   d_directory_name = directory_name;
}

/*
*************************************************************************
*                                                                       *
* Register standard (non-derived) plot quantities, vector or scalar.    *
*                                                                       *
*************************************************************************
*/
template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerPlotScalar(
   const std::string& variable_name, 
   int data_id, 
   int depth_id,
   double scale_factor)
{
   registerPlotItem(variable_name, data_id, false, depth_id, scale_factor);
}

template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerPlotVector(
   const std::string& variable_name,
   int data_id,
   double scale_factor)
{
   registerPlotItem(variable_name, data_id, true, 0, scale_factor);
}

/*
*************************************************************************
*                                                                       *
* Private function to register a standard (non-derived) plot quantity.  *
* We check to make sure that the variable name is unique, that the      *
* factory at given index is defined and that the data is cell-centered. *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerPlotItem(
   const std::string& variable_name, 
   int data_id,
   bool isa_vector,
   int depth_id,
   double scale_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(data_id > -1);
   TBOX_ASSERT(depth_id > -1);
#endif

   for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items);
        ipi; ipi++) {
      if (ipi().d_variable_name == variable_name) {
         TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::registerPlotItem"
         << " error...\n"
         << "    Attempting to register variable with name " << variable_name
         << "\n     more than once." << std::endl);
      }
   }

   VizamraiItem<DIM> plotitem;
   plotitem.d_variable_name   = variable_name; 
   plotitem.d_data_id         = data_id; 
   plotitem.d_isa_vector      = isa_vector;
   plotitem.d_depth           = 1; // default to 1, if vector lookup below
   plotitem.d_depth_id        = depth_id; 
   plotitem.d_scale_factor    = scale_factor; 
   plotitem.d_derived_writer  = NULL; 

   tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
      hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
                                             getPatchDataFactory(data_id);
   if (factory.isNull()) {
      TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::registerPlotItem"
         << " error...\n"
         << "    data id = " << data_id << " for variable = " << variable_name
         << "\n     is undefined in the patch descriptor" << std::endl);
   } else {
      bool found_type = false;
      bool depth_is_good = false;
#ifdef HAVE_FLOAT
      tbox::Pointer< pdat::CellDataFactory<DIM,float> > ffactory = factory;
      if (!ffactory.isNull()) { 
         found_type = true;
	 if (isa_vector) {
	    plotitem.d_depth  = ffactory->getDefaultDepth();
            depth_is_good = true;
         } else {
            depth_is_good = ( depth_id < ffactory->getDefaultDepth() );
         }
         plotitem.d_data_type = DATA_TO_VIZ_IS_FLOAT;
      }
#endif
      if (!found_type) {
         tbox::Pointer< pdat::CellDataFactory<DIM,double> > dfactory = factory;
         if (!dfactory.isNull()) { 
            found_type = true;
	    if (isa_vector) {
	       plotitem.d_depth  = dfactory->getDefaultDepth();
               depth_is_good = true;
            } else {
               depth_is_good = ( depth_id < dfactory->getDefaultDepth() );
            }
            plotitem.d_data_type = DATA_TO_VIZ_IS_DOUBLE;      
         }
      }
      if (!found_type) {
         tbox::Pointer< pdat::CellDataFactory<DIM,int> > ifactory = factory;
         if (!ifactory.isNull()) { 
            found_type = true;
	    if (isa_vector) {
	       plotitem.d_depth  = ifactory->getDefaultDepth();
               depth_is_good = true;
            } else {
               depth_is_good = ( depth_id < ifactory->getDefaultDepth() );
            }
            plotitem.d_data_type = DATA_TO_VIZ_IS_INT;
         }
      }
      if (!found_type) {
         TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::registerPlotItem"
           << " error...\n"
           << "    data id = " << data_id << " for variable = " 
           << variable_name 
           << "\n     corresponds to an illegal patch descriptor entry." 
           << "\n     Legal entries are cell-centered double, float, or int."
           << std::endl);
      }
      if (!depth_is_good) {
         TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::registerPlotItem"
           << " error...\n"
           << "    depth id = " << depth_id << " for variable = "
           << variable_name
           << "\n     is out of range." << std::endl);
      }
   }

   d_plot_items.appendItem(plotitem);
}

/*
*************************************************************************
*                                                                       *
* Reset previously-registered non-derived plot quantity to new data id  *
* and depth id on the given level.  This allows us to use different     * 
* patch data ids for the same quantity on different hierarchy levels.   *
* We check to make sure that the factory at the given index is          *
* defined and consistent with the original registration.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::resetLevelPlotScalar(
   const std::string& variable_name,
   int level_number,
   int data_id,
   int depth_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(level_number > -1);
   TBOX_ASSERT(data_id > -1);
   TBOX_ASSERT(depth_id > -1);
#endif

   bool found = false;
   for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items);
        (!found && ipi); ipi++) {
      if (ipi().d_variable_name == variable_name) {

         tbox::Pointer< hier::PatchDataFactory<DIM> > factory =
         hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
                                                getPatchDataFactory(data_id);
         

         if (factory.isNull()) {
            TBOX_ERROR(
               "CartesianVizamraiDataWriter<DIM>::resetLevelPlotScalar"
               << " error...\n"
               << "    data id = " << data_id << " for variable = " 
               << variable_name
               << "\n     is undefined in the patch descriptor" << std::endl);
         }
    
         int old_data_type = ipi().d_data_type;

         bool depth_is_good = false;
         bool found_type = false;

#ifdef HAVE_FLOAT 
         if (old_data_type == DATA_TO_VIZ_IS_FLOAT) {
            tbox::Pointer< pdat::CellDataFactory<DIM,float> > ffactory = factory;
            if (!ffactory.isNull()) {
               found_type = true;
               depth_is_good = ( depth_id < ffactory->getDefaultDepth() );
            }
         } 
#endif

         if (old_data_type == DATA_TO_VIZ_IS_DOUBLE) {
            tbox::Pointer< pdat::CellDataFactory<DIM,double> > dfactory = factory;
            if (!dfactory.isNull()) {
               found_type = true;
               depth_is_good = ( depth_id < dfactory->getDefaultDepth() );
            }
         }

         if (old_data_type == DATA_TO_VIZ_IS_INT) {
            tbox::Pointer< pdat::CellDataFactory<DIM,int> > ifactory = factory;
            if (!ifactory.isNull()) {
               found_type = true;
               depth_is_good = ( depth_id < ifactory->getDefaultDepth() );
            }
         }

         if (!found_type) {
            TBOX_ERROR(
              "CartesianVizamraiDataWriter<DIM>::resetLevelPlotScalar"
              << " error...\n"
              << "    data id = " << data_id << " for variable = "
              << variable_name
              << "\n     does not match patch data type of previously" 
              << "\n     registered variable" << std::endl);
         }
         if (!depth_is_good) {
            TBOX_ERROR(
              "CartesianVizamraiDataWriter<DIM>::resetLevelPlotScalar"
              << " error...\n"
              << "    depth id = " << depth_id << " for variable = "
              << variable_name << "\n     is out of range." << std::endl);
         }

         found = true; 

         tbox::Array<int>& plot_data_id = ipi().d_level_data_id;
         tbox::Array<int>& plot_depth_id = ipi().d_level_depth_id;

         int oldsize = plot_data_id.getSize();
         int newsize = level_number + 1;
         if ( oldsize < newsize ) {
            plot_data_id.resizeArray(newsize);
            plot_depth_id.resizeArray(newsize);
            for (int i = oldsize; i < newsize; i++) {
               plot_data_id[i] = DATA_TO_VIZ_IS_UNDEFINED;
               plot_depth_id[i] = DATA_TO_VIZ_IS_UNDEFINED;
            }
         }

         plot_data_id[level_number] = data_id;
         plot_depth_id[level_number] = depth_id;

      } // test on variable name
   } // iteration over list of plot items

   if (!found) {
      TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::resetLevelPlotScalar"
                 << " error...\n"
                 << "    variable = " << variable_name
                 << " has not been registered for plotting." << std::endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Reset previously-registered non-derived plot quantity to new data id  *
* and depth id on the given level.  This allows us to use different     * 
* patch data ids for the same quantity on different hierarchy levels.   *
* We check to make sure that the factory at the given index is          *
* defined and consistent with the original registration.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::resetLevelPlotVector(
   const std::string& variable_name,
   int level_number,
   int data_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(level_number > -1);
   TBOX_ASSERT(data_id > -1);
#endif

   bool found = false;
   for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items);
        (!found && ipi); ipi++) {
      if (ipi().d_variable_name == variable_name) {

         tbox::Pointer< hier::PatchDataFactory<DIM> > factory =
         hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
                                                getPatchDataFactory(data_id);
         
         if (factory.isNull()) {
            TBOX_ERROR(
               "CartesianVizamraiDataWriter<DIM>::resetLevelPlotVector"
               << " error...\n"
               << "    data id = " << data_id << " for variable = " 
               << variable_name
               << "\n     is undefined in the patch descriptor" << std::endl);
         }
    

         int old_data_type = ipi().d_data_type;
         int old_data_depth = ipi().d_depth;
         int new_data_depth = 0;
         
         bool depth_is_good = false;
         bool found_type = false;

#ifdef HAVE_FLOAT 
         if (old_data_type == DATA_TO_VIZ_IS_FLOAT) {
            tbox::Pointer< pdat::CellDataFactory<DIM,float> > ffactory = factory;
            if (!ffactory.isNull()) {
               found_type = true;
               new_data_depth = ffactory->getDefaultDepth();
               depth_is_good = ( old_data_depth == new_data_depth );
            }
         } 
#endif

         if (old_data_type == DATA_TO_VIZ_IS_DOUBLE) {
            tbox::Pointer< pdat::CellDataFactory<DIM,double> > dfactory = factory;
            if (!dfactory.isNull()) {
               found_type = true;
               new_data_depth = dfactory->getDefaultDepth();
               depth_is_good = ( old_data_depth == new_data_depth );
            }
         }

         if (old_data_type == DATA_TO_VIZ_IS_INT) {
            tbox::Pointer< pdat::CellDataFactory<DIM,int> > ifactory = factory;
            if (!ifactory.isNull()) {
               found_type = true;
               new_data_depth = ifactory->getDefaultDepth();
               depth_is_good = ( old_data_depth == new_data_depth );
            }
         }

         if (!found_type) {
            TBOX_ERROR(
              "CartesianVizamraiDataWriter<DIM>::resetLevelPlotVector"
              << " error...\n"
              << "    data id = " << data_id << " for variable = "
              << variable_name
              << "\n     does not match patch data type of previously" 
              << "\n     registered variable" << std::endl);
         }
         if (!depth_is_good) {
            TBOX_ERROR(
              "CartesianVizamraiDataWriter<DIM>::resetLevelPlotVector"
              << " error...\n"
              << "    depth of new variable = " << new_data_depth
              << "\n    depth of old variable = " << old_data_depth
              << std::endl);
         }

         found = true; 

         tbox::Array<int>& plot_data_id = ipi().d_level_data_id;

         int oldsize = plot_data_id.getSize();
         int newsize = level_number + 1;
         if ( oldsize < newsize ) {
            plot_data_id.resizeArray(newsize);
            for (int i = oldsize; i < newsize; i++) {
               plot_data_id[i] = DATA_TO_VIZ_IS_UNDEFINED;
            }
         }

         plot_data_id[level_number] = data_id;

      } // test on variable name
   } // iteration over list of plot items

   if (!found) {
      TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::resetLevelPlotVector"
                 << " error...\n"
                 << "    variable = " << variable_name
                 << " has not been registered for plotting." << std::endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Register derived plot quantities, scalar or vector.  If no derived    *
* data strategy is specified and no default derived data strategy is    * 
* set, an error will result.                                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerDerivedPlotScalar(
   const std::string& variable_name,
   VisDerivedDataStrategy<DIM>* derived_writer) 
{
   registerDerivedPlotItem(variable_name, 1, derived_writer);
}

template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerDerivedPlotVector(
   const std::string& variable_name,
   int depth,
   VisDerivedDataStrategy<DIM>* derived_writer)
{
   registerDerivedPlotItem(variable_name, depth, derived_writer);
}

/*
*************************************************************************
*                                                                       *
* Private function to register a derived plot quantity.  We check to    *
* make sure that the variable name is unique and that a derived data    *
* strategy object is available.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::registerDerivedPlotItem(
   const std::string& variable_name,
   int depth,
   VisDerivedDataStrategy<DIM>* derived_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(depth > 0);
#endif

   for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items); 
        ipi; ipi++) {
      if (ipi().d_variable_name == variable_name) {
         TBOX_ERROR("CartesianVizamraiDataWriter<DIM>::registerDerivedPlotItem"
         << " error...\n"
         << "    Attempting to register variable with name " << variable_name
         << "\n     more than once." << std::endl);
      }
   }

   bool use_default_writer = (derived_writer == NULL);

   if (use_default_writer && (d_default_derived_writer == NULL)) {
      TBOX_ERROR(
         "CartesianVizamraiDataWriter<DIM>::registerDerivedPlotItem"
         << " error...\n"
         << "    no derived data writer specified for variable = "
         << variable_name
         << "\n    and no default derived data writer set." << std::endl);
   } 

   VizamraiItem<DIM> plotitem;
   plotitem.d_variable_name   = variable_name;
   plotitem.d_data_id         = DATA_TO_VIZ_IS_UNDEFINED;
   plotitem.d_isa_vector      = (depth > 1 ? true : false);
   plotitem.d_depth_id        = DATA_TO_VIZ_IS_UNDEFINED;
   plotitem.d_depth           = depth;
   plotitem.d_scale_factor    = 1.0;
   plotitem.d_derived_writer  = (use_default_writer ? d_default_derived_writer
                                                    : derived_writer);
   plotitem.d_data_type       = DATA_TO_VIZ_IS_UNDEFINED;

   d_plot_items.appendItem(plotitem);
}

/*
*************************************************************************
*                                                                       *
* Write plot data from given hierarchy to file with given name.         * 
* If no directory is specified, the default directory is used.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::writePlotData(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const std::string& file_name,
   int extension,
   double plot_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT(!file_name.empty());
#endif 

   /*
    * First create a temporary file, when all data has been written rename it.
    */
   tbox::FileStream *file = 
      new tbox::FileStream(
          (getFileStreamName(file_name, extension, true)).c_str(), 
           tbox::FileStream::Write);

   /*
    * Determine finest level to plot, adjusting parameters if appropriate.
    */ 

   int finest_level_to_plot = 
      tbox::MathUtilities<int>::Min(hierarchy->getFinestLevelNumber(),
                                    d_finest_plot_level_now);

   if (!d_wrote_plot_file && (d_finest_plot_level_in_all < 0)) {
      finest_level_to_plot = hierarchy->getFinestLevelNumber(); 
      d_finest_plot_level_in_all = finest_level_to_plot;
      d_finest_plot_level_now = finest_level_to_plot;
   } 

   /*
    * Check whether all plot levels exists and if level scaling information
    * matches ratios in AMR hierarchy.  If everything is good, we write
    * plot file; otherwise, we exit.
    */ 

   if (checkLevelInformation(hierarchy, 0, finest_level_to_plot)) {

      /*
       * Compute box regions on each level to write to plot file; 
       * removing refined regions as needed.
       */

      tbox::Array< hier::BoxList<DIM> >* output_boxes = new 
         tbox::Array< hier::BoxList<DIM> >[finest_level_to_plot+1];
      int num_outboxes = computeOutputBoxes(output_boxes, 
                                            hierarchy,
                                            0, finest_level_to_plot);

      /*
       * Compute scaled domain length used when DIM < 3.
       */

      const tbox::Pointer< geom::CartesianGridGeometry<DIM> > ggeom = 
         hierarchy->getGridGeometry();
      const double* const xupper = ggeom->getXUpper();
      const double* const xlower = ggeom->getXLower();
 
      d_domain_scale_length = (xupper[0] - xlower[0]);
      for (int id = 1; id < DIM; id++) {
         d_domain_scale_length = 
            tbox::MathUtilities<double>::Min( d_domain_scale_length,
                                              (xupper[id] - xlower[id]) );
      }
      d_domain_scale_length *= DOMAIN_SCALE_FACTOR;

      /*
       * Write Vizamrai header file and patch boundary information.
       */

      writeVizamraiHeaderInfoToFile(*file,
                                    plot_time,
                                    num_outboxes);

      writePatchBoundariesToFile(*file,
                                 hierarchy,
                                 0, finest_level_to_plot);

      /*
       * Write visualization data to plot file.
       */

      writeVizamraiVariablesToFile(*file,
                                   hierarchy,
                                   0, finest_level_to_plot,
                                   output_boxes); 

      delete[] output_boxes;

      d_wrote_plot_file = true;

      delete file;

      /*
       * Move the temporary file to the final destination.
       */
      if (tbox::SAMRAI_MPI::getNodes() == 1) {
	tbox::Utilities::renameFile(
           (getFileStreamName(file_name, extension, true)).c_str(),
	   (getFileStreamName(file_name, extension, false)).c_str());
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Print all class data for HyperbolicLevelIntegrator object.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::printClassData(std::ostream& os) const
{
   os << "\nCartesianVizamraiDataWriter<DIM>::printClassData..." << std::endl;
   os << "CartesianVizamraiDataWriter<DIM>: this = "
      << (CartesianVizamraiDataWriter<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_wrote_plot_file = " << d_wrote_plot_file << std::endl;
   os << "d_finest_plot_level_in_all = " << d_finest_plot_level_in_all << std::endl;
   os << "d_finest_plot_level_now = " << d_finest_plot_level_now << std::endl;
   for (int i = 0; i < d_scaling_ratios.getSize(); i++) {
      os << "d_scaling_ratios[" << i << "] = " 
         <<  d_scaling_ratios[i] << std::endl;
   }
   os << "d_plot_type = "; 
   if (d_plot_type == DUMP_FLOAT_DATA) {
      os << "Float Data" << std::endl;
   } else {
      os << "Double Data" << std::endl;
   } 
   os << "d_default_derived_writer = " << d_default_derived_writer << std::endl;
   os << "d_directory_name = " << d_directory_name << std::endl;
   os << "d_plot_items list..." << std::endl; 
   os << "Number of plot items = " 
      << d_plot_items.getNumberItems() << std::endl; 
   for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items); 
        ipi; ipi++) {
      os << "   Name = " << ipi().d_variable_name << std::endl;
      if (!(ipi().d_derived_writer == NULL)) {
         if (ipi().d_isa_vector) {
            os << "     DERIVED PLOT VECTOR: depth = " << ipi().d_depth << std::endl;   
         } else {
            os << "     DERIVED PLOT SCALAR" << std::endl;   
         }
      } else {
         os << "     Master hier::Patch Data ID = " << ipi().d_data_id << std::endl;
         for (int kk = 0; kk < ipi().d_level_data_id.getSize(); kk++) {
            os << "        Data id on level " << kk << " is " 
                                                    << ipi().d_level_data_id[kk] << std::endl;
         }
         if (ipi().d_isa_vector) {
            os << "     Data is a vector" << std::endl; 
            os << "        Depth is" << ipi().d_depth << std::endl; 
         } else {
            os << "     Data is a scalar" << std::endl;    
            os << "        Master hier::Patch Data Depth ID = " << ipi().d_depth_id << std::endl;
            for (int ll = 0; ll < ipi().d_level_data_id.getSize(); ll++) {
               os << "        Depth id on level " << ll << " is "
                                                        << ipi().d_level_depth_id[ll] << std::endl;
            }
         }
         os << "     Scale Factor = " << ipi().d_scale_factor << std::endl;
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Private utility function to parse file name, directory name, and      *
* extension and return an appropriate file stream name.  Note that when *
* running on a single processor, we don't split output files.  When     *
* running in parallel, node 0 outputs global information and output     *
* from each node is written to a seperate file that must later be       *
* concatenated for viewing.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> std::string CartesianVizamraiDataWriter<DIM>::getFileStreamName(
   const std::string& file_name,
   int extension,
   bool istemporaryflag) const
{
   std::string dump_directory = d_directory_name;

   std::string dump_filename = file_name;
   if (!dump_directory.empty()) {
      tbox::Utilities::recursiveMkdir(dump_directory);
      std::string tmp = dump_filename;
      dump_filename = dump_directory;
      dump_filename += "/"; 
      dump_filename += tmp;
   }

   const int size = dump_filename.length() + 20;
   char *buffer = new char[size];

   if (tbox::SAMRAI_MPI::getNodes() > 1) {
      if (extension >= 0) {
         sprintf(buffer, "%s.%05d.vis.%05d", dump_filename.c_str(),
                 extension, tbox::SAMRAI_MPI::getRank());
      } else {
         sprintf(buffer, "%s.vis.%05d", dump_filename.c_str(),
                 tbox::SAMRAI_MPI::getRank());
      }
   } else {
      // If this is a temporary file then add an tmp extension to the name
      if(istemporaryflag) {
	 if (extension >= 0) {
	    sprintf(buffer, "%s.%05d.vis.tmp", dump_filename.c_str(), extension);
	 } else {
	    sprintf(buffer, "%s.vis.tmp", dump_filename.c_str());
	 }
      } else {
	 if (extension >= 0) {
	    sprintf(buffer, "%s.%05d.vis", dump_filename.c_str(), extension);
	 } else {
	    sprintf(buffer, "%s.vis", dump_filename.c_str());
	 }
      }
   }

   std::string stream_name(buffer); 

   delete [] buffer;

   return(stream_name);
}

/*
*************************************************************************
*                                                                       *
* Private utility function that writes Vizamrai header information.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::writeVizamraiHeaderInfoToFile(
   tbox::FileStream& file,
   double plot_time,
   int num_outboxes) const
{
   if (tbox::SAMRAI_MPI::getRank() == 0) {

      /*
       * Vizamrai header information:
       *    - Vizamrai file format version
       *    - timestamp for plot data
       *    - number of box regions in data file
       *    - type of plot data (float or double)
       *    - number of plot quantities
       *    - string identifier for each plot quantity
       *    - maximum number of levels over all plot files
       *    - mesh ratio between each level and coarsest level (i.e., level 0)
       */

      file << -VIZAMRAI_FILE_FORMAT_VERSION;  
      file << plot_time;
      file << num_outboxes;
      file << d_plot_type;

      file << d_plot_items.getNumberItems();

      for (typename tbox::List<VizamraiItem<DIM> >::Iterator ipi(d_plot_items); 
           ipi; ipi++) {
         file.writeString(ipi().d_variable_name.c_str());
	 file << ipi().d_depth;
      }

      file << (d_finest_plot_level_in_all + 1); 

      for (int ln = 1; ln <= d_finest_plot_level_in_all; ln++) {
         const int* tmpvec = d_scaling_ratios[ln];
         file.pack(tmpvec, DIM);
	 if (DIM < 3) {
	    file << 1;
	 }
	 if (DIM == 1) {
	    file << 1;
	 }
      } 
      
   }
}

/*
*************************************************************************
*                                                                       *
* Private utility function to compute box regions on each level         *
* containing data that will be written to Vizamrai plot file.           *
* That is, this routine removes from each level regions that are        *
* covered by some finer level.  We return the total number of box       *
* regions containing plot data.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> int CartesianVizamraiDataWriter<DIM>::computeOutputBoxes(
   tbox::Array< hier::BoxList<DIM> > outboxes[],
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int coarsest_plot_level,
   int finest_plot_level) const
{
   int num_outboxes = 0;

   for (int ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {

      tbox::Pointer<hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
      const hier::BoxArray<DIM>& level_boxes = level->getBoxes();
      const int num_patches = level_boxes.getNumberOfBoxes();

      hier::BoxList<DIM> tmp_outboxes = level_boxes;
      if (ln < finest_plot_level) {
         tbox::Pointer<hier::PatchLevel<DIM> > finer_level =
            hierarchy->getPatchLevel(ln+1);
         hier::IntVector<DIM> coarsen_ratio = finer_level->getRatioToCoarserLevel();
         hier::BoxList<DIM> takeaway = finer_level->getBoxes();
         takeaway.coarsen(coarsen_ratio);
         tmp_outboxes.removeIntersections(takeaway);
      }

      outboxes[ln].resizeArray(num_patches);

      for (int ilb = 0; ilb < num_patches; ilb++) {
         hier::BoxList<DIM>& patch_outboxes = outboxes[ln][ilb];

         for (typename hier::BoxList<DIM>::Iterator iob(tmp_outboxes); iob; iob++) {
            hier::Box<DIM> intersection = level_boxes(ilb) * iob(); 
            if (!intersection.empty()) {
               patch_outboxes.appendItem(intersection);
               num_outboxes++;
            }
         }

      }
             
   }

   return(num_outboxes);
}

/*
*************************************************************************
*                                                                       *
* Private utility function that writes patch box information to file.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::writePatchBoundariesToFile(
   tbox::FileStream& file,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int coarsest_plot_level,
   int finest_plot_level) const
{
   int ln, id;

   int noutput = 0;

   if (tbox::SAMRAI_MPI::getRank() == 0) {

      for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
         tbox::Pointer<hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
         const hier::BoxArray<DIM>& boxes = level->getBoxes();
         noutput += boxes.getNumberOfBoxes();
      }
      file << noutput;

   }

   for (ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         const tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());
         const tbox::Pointer< geom::CartesianPatchGeometry<DIM> > geometry =
            patch->getPatchGeometry();
         const double* xlo = geometry->getXLower();
         const double* xup = geometry->getXUpper();
 
         if (tbox::SAMRAI_MPI::getRank() == 0) {
 
            noutput--;
 
            for (id = 0; id < DIM; id++) {
               file << xlo[id];
            }
	    if (DIM < 3) {
	       file  << 0.0; 
	    }
	    if (DIM==1) {
	       file  << 0.0;
	    }

            for (id = 0; id < DIM; id++) {
               file << xup[id];
            }
	    if (DIM < 3) {
	       file  << d_domain_scale_length;
	    }
	    if (DIM==1) {
	       file  << d_domain_scale_length;
	    }
            file << ln;
 
         } else {

            double buffer[DIM*2+1];
            int buffer_count = 0;

            for (id = 0; id < DIM; id++) {
               buffer[buffer_count++] = xlo[id];
            }

            for (id = 0; id < DIM; id++) {
               buffer[buffer_count++] = xup[id];
            }

#ifdef HAVE_MPI
            buffer[buffer_count++] = ln;
            MPI_Send(buffer, DIM*2+1, MPI_DOUBLE,
                     0, 0, tbox::SAMRAI_MPI::getCommunicator());
#endif
         }

      } // loop over patches 

   } // loop over levels

#ifdef HAVE_MPI
   MPI_Status status;
#endif

   if (tbox::SAMRAI_MPI::getRank() == 0) {

      while(noutput) {
   
         double buffer[DIM*2+1];

         noutput--;

#ifdef HAVE_MPI
         MPI_Recv(buffer, DIM*2+1, MPI_DOUBLE, MPI_ANY_SOURCE, 0,
                  tbox::SAMRAI_MPI::getCommunicator(), &status);
#endif

         for (id = 0; id < DIM; id++) {
            file << buffer[id];
         }
	 if (DIM < 3) {
	    file  << 0.0;
	 }
	 if (DIM == 1) {
	    file  << 0.0; 
	 }

         for (id = 0; id < DIM; id++) {
            file << buffer[id+DIM];
         }

	 if (DIM < 3) {
	    file  << d_domain_scale_length;
	 }
	 if (DIM == 1) {
	    file  << d_domain_scale_length;
	 }

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
         file << static_cast<int>(buffer[DIM*2]);

      }

   } 
}

/*
*************************************************************************
*                                                                       *
* Private utility function to check scaling ratio data between levels.  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::writeVizamraiVariablesToFile(
   tbox::FileStream& file,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int coarsest_plot_level,
   int finest_plot_level,
   const tbox::Array< hier::BoxList<DIM> > outboxes[]) const
{

   const tbox::Pointer< geom::CartesianGridGeometry<DIM> > ggeom =
      hierarchy->getGridGeometry();
   const double* xlo_domain = ggeom->getXLower();
   const double* dx_domain = ggeom->getDx();

   tbox::DatabaseBox bounding_box = hier::BoxList<DIM>(ggeom -> getPhysicalDomain()).getBoundingBox();

   // Center of the 0,0,0 cell
   double xzero_domain[DIM];
   for (int id = 0; id < DIM; id++) {
      xzero_domain[id] = xlo_domain[id] - bounding_box.lower(id) * dx_domain[id];
   }

   for (int ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);

      for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip()); 

         const tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pgeom =
            patch->getPatchGeometry();
         const double* const dx_patch = pgeom->getDx();
	 
         const hier::BoxList<DIM>& patch_plot_boxes = outboxes[ln][ip()];

         for (typename hier::BoxList<DIM>::Iterator ipb(patch_plot_boxes); ipb; ipb++) {

            writePlotBoxDataToFile(file,
                                   ln,
                                   xzero_domain,
                                   dx_patch,
                                   ipb());
 
            const int bsize = ipb().size();

            double* dbuffer = new double[bsize];
            float*  fbuffer = NULL; 
            if (d_plot_type == DUMP_FLOAT_DATA) {
               fbuffer = new float[bsize];
            }

            for (typename tbox::List<VizamraiItem<DIM> >::Iterator 
                 ipi(d_plot_items); ipi; ipi++) {
               
               if (!(ipi().d_derived_writer == NULL)) {

                  /*
                   * Write derived patch data
                   */
                  int depth_id = ipi().d_depth_id;
                  if ((ln < ipi().d_level_data_id.getSize()) &&
                      (ipi().d_level_data_id[ln] > DATA_TO_VIZ_IS_UNDEFINED)) {
                     depth_id = ipi().d_level_depth_id[ln];
                  }

		  if (ipi().d_isa_vector) {

		     int depth = ipi().d_depth;
		     for(depth_id = 0; depth_id < depth; depth_id ++) {

                        (void) ipi().d_derived_writer->
                           packDerivedDataIntoDoubleBuffer(
                              dbuffer,
                              *patch,
                              ipb(),
                              ipi().d_variable_name, 
                              depth_id);

                        if (d_plot_type == DUMP_FLOAT_DATA) {
                           for (int i = 0; i < bsize; i++) {
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
                              fbuffer[i] = static_cast<float>(dbuffer[i]);
                           }
                           file.pack(fbuffer, bsize);
                        } else {
                           file.pack(dbuffer, bsize);
                        }
                     }
                     
                  } else {
                     
                     (void) ipi().d_derived_writer->
                        packDerivedDataIntoDoubleBuffer(
                           dbuffer,
                           *patch,
                           ipb(),
                           ipi().d_variable_name, 
                           depth_id);

                     if (d_plot_type == DUMP_FLOAT_DATA) {
                        for (int i = 0; i < bsize; i++) {
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
                           fbuffer[i] = static_cast<float>(dbuffer[i]);
                        }
                        file.pack(fbuffer, bsize);
                     } else {
                        file.pack(dbuffer, bsize);
                     }
                  }

               } else {

                  int data_id = ipi().d_data_id;
                  int depth_id = ipi().d_depth_id;
                  if ((ln < ipi().d_level_data_id.getSize()) &&
                      (ipi().d_level_data_id[ln] > DATA_TO_VIZ_IS_UNDEFINED)) {
                     data_id = ipi().d_level_data_id[ln];
                     depth_id = ipi().d_level_depth_id[ln];
                  }

		  if (ipi().d_isa_vector) {

		     int depth = ipi().d_depth;
		     for(depth_id = 0; depth_id < depth; depth_id ++) {

			packPatchDataIntoDoubleBuffer(
			   patch->getPatchData(data_id),
			   depth_id,
			   ipi().d_data_type,
			   ipb(),
			   dbuffer);                  
		     
			const double scale = ipi().d_scale_factor;
// Disable Intel warning about real comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
			if (scale != 1.0) {
			   for (int i = 0; i < bsize; i++) {
			      dbuffer[i] *= scale;
			   }
			}
			
			if (d_plot_type == DUMP_FLOAT_DATA) {
			   for (int i = 0; i < bsize; i++) {
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
			      fbuffer[i] = static_cast<float>(dbuffer[i]);
			   }
			   file.pack(fbuffer, bsize);
			} else {
			   file.pack(dbuffer, bsize);
			}

		     }

		  } else {

		     packPatchDataIntoDoubleBuffer(
			patch->getPatchData(data_id),
			depth_id,
			ipi().d_data_type,
			ipb(),
			dbuffer);                  
		     
		     const double scale = ipi().d_scale_factor;
// Disable Intel warning about real comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
		     if (scale != 1.0) {
			for (int i = 0; i < bsize; i++) {
			   dbuffer[i] *= scale;
			}
		     }
		     
		     if (d_plot_type == DUMP_FLOAT_DATA) {
			for (int i = 0; i < bsize; i++) {
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
			   fbuffer[i] = static_cast<float>(dbuffer[i]);
			}
			file.pack(fbuffer, bsize);
		     } else {
			file.pack(dbuffer, bsize);
		     }
		  }

               }  // if not derived plot data

            } // loop over plot items

            if (d_plot_type == DUMP_FLOAT_DATA) {
               delete [] fbuffer;
            }
            delete [] dbuffer;

         }  // loop over plot boxes on patch

      }  // loop over patches

   }  // loop over levels
}

/*
*************************************************************************
*                                                                       *
* Private utility function to write single box information to file.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::writePlotBoxDataToFile(
   tbox::FileStream& file,
   int level_number,
   const double* const xzero_domain,
   const double* const dx_level,
   const hier::Box<DIM>& plot_box) const
{
   file << level_number;

   for (int id1 = 0; id1 < DIM; id1++) {
      file << plot_box.lower(id1);
   }
   if (DIM < 3) {
      file << 0;
   }
   if (DIM == 1) {
      file << 0;
   }

   for (int id2 = 0; id2 < DIM; id2++) {
      file << plot_box.upper(id2);
   }
   if (DIM < 3) {
      file << 0;
   }
   if (DIM == 1) {
      file << 0;
   }

   for (int id3 = 0; id3 < DIM; id3++) {
      file << dx_level[id3];
   }
   if (DIM < 3) {
      file << d_domain_scale_length;
   }
   if (DIM == 1) {
      file << d_domain_scale_length;
   }

   for (int id4 = 0; id4 < DIM; id4++) {
      const double level_domain_zero_cell_center = 
         xzero_domain[id4] + 0.5 * dx_level[id4];
      file << level_domain_zero_cell_center;
   }
   if (DIM < 3) {
      file << (d_domain_scale_length/2.0);
   }
   if (DIM == 1) {
      file << (d_domain_scale_length/2.0);
   }

}

/*
*************************************************************************
*                                                                       *
* Private utility function to check scaling ratio data between levels.  *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool CartesianVizamraiDataWriter<DIM>::checkLevelInformation(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int coarsest_plot_level,
   int finest_plot_level)
{
   bool ret_val = true;

   for (int ln = coarsest_plot_level; ln <= finest_plot_level; ln++) {
      if (hierarchy->getPatchLevel(ln).isNull()) {
         ret_val = false;
         tbox::perr << "CartesianVizamraiDataWriter<DIM>::writePlotData error...\n";
         tbox::perr << "ERROR MESSAGE: patch level " << ln << " is null." << std::endl;
         printClassData(tbox::perr);
         tbox::perr << std::flush;
         tbox::SAMRAI_MPI::abort();
      }
   }

   int coarsest_to_check = 
      tbox::MathUtilities<int>::Max(1, coarsest_plot_level);

   if ( !d_wrote_plot_file && (d_scaling_ratios.getSize() == 0) ) {
      for (int ln = coarsest_to_check; ln <= finest_plot_level; ln++) {
         tbox::Pointer<hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
         d_scaling_ratios[ln] = level->getRatioToCoarserLevel();
      }
   } else {
      if (finest_plot_level+1 > d_scaling_ratios.getSize()) {
         ret_val = false; 
      } else {
         for (int ln = coarsest_to_check; 
              ret_val && (ln <= finest_plot_level); ln++) {
            tbox::Pointer<hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
            if ( !(d_scaling_ratios[ln] == level->getRatioToCoarserLevel()) ) {
               ret_val = false;
            }
         }
      }
   }

   if (ret_val == false) {
      tbox::perr << "CartesianVizamraiDataWriter<DIM>::writePlotData error...\n";
      tbox::perr << "ERROR MESSAGE: Level ratios set in setRatioToCoarserLevel()";
      tbox::perr << "\n do not match ratios between levels in hierarchy." << std::endl;
      printClassData(tbox::perr);
      tbox::perr << std::flush;
      tbox::SAMRAI_MPI::abort();
   }

   return(ret_val);
}

/*
*************************************************************************
*                                                                       *
* Private utility function to pack data into double buffer.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CartesianVizamraiDataWriter<DIM>::packPatchDataIntoDoubleBuffer(
   const tbox::Pointer< hier::PatchData<DIM> > pdata,
   int depth,
   int type_of_data,
   const hier::Box<DIM>& box,
   double* buffer) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((box * pdata->getGhostBox()) == box);
#endif

   const hier::Box<DIM>& data_box = pdata->getGhostBox();

   const int box_w0 = box.numberCells(0);
   int dat_w0 = -1;
   int box_w1 = -1;

   if (DIM > 1) {
      dat_w0 = data_box.numberCells(0);
      box_w1 = box.numberCells(1);
   }

   int dat_w1 = -1;
   int box_w2 = -1;

   if (DIM > 2) {
      dat_w1 = data_box.numberCells(1);
      box_w2 = box.numberCells(2);
   }

   int buf_b1 = 0;

   int dat_b2 = data_box.offset(box.lower());

   double *const buf_ptr = buffer;

   switch(type_of_data) {

      case DATA_TO_VIZ_IS_FLOAT: {
#ifdef HAVE_FLOAT
         const tbox::Pointer< pdat::CellData<DIM,float> > fpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!fpdata.isNull());
#endif
         const float *const dat_ptr = fpdata->getPointer(depth);


	 if (DIM == 1) {
            int dat_b1 = dat_b2;
	    for (int i0 = 0; i0 < box_w0; i0++) {
	       buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
	    }
	 } else if (DIM == 2) {
            int dat_b1 = dat_b2;
            for (int i1 = 0; i1 < box_w1; i1++) {
               for (int i0 = 0; i0 < box_w0; i0++) {
                  buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
               }
               dat_b1 += dat_w0;
               buf_b1 += box_w0;
            }
	 } else if (DIM == 3) {
	    for (int i2 = 0; i2 < box_w2; i2++) {
	       int dat_b1 = dat_b2;
	       for (int i1 = 0; i1 < box_w1; i1++) {
		  for (int i0 = 0; i0 < box_w0; i0++) {
		     buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
		  }
		  dat_b1 += dat_w0;
		  buf_b1 += box_w0;
	       }
	       dat_b2 += dat_w1 * dat_w0;
	    }
	 }

#endif
         break;

      }

      case DATA_TO_VIZ_IS_DOUBLE: {

         const tbox::Pointer< pdat::CellData<DIM,double> > dpdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!dpdata.isNull());
#endif
         const double *const dat_ptr = dpdata->getPointer(depth);

	 if (DIM == 1) {
	    int dat_b1 = dat_b2;
	    for (int i0 = 0; i0 < box_w0; i0++) {
	       buf_ptr[buf_b1+i0] = dat_ptr[dat_b1+i0];
	    }
	 } else if (DIM == 2) {
            int dat_b1 = dat_b2;
            for (int i1 = 0; i1 < box_w1; i1++) {
               for (int i0 = 0; i0 < box_w0; i0++) {
                  buf_ptr[buf_b1+i0] = dat_ptr[dat_b1+i0];
               }
               dat_b1 += dat_w0;
               buf_b1 += box_w0;
            }
	 } else if (DIM == 3) {
	    for (int i2 = 0; i2 < box_w2; i2++) {
	       int dat_b1 = dat_b2;
	       for (int i1 = 0; i1 < box_w1; i1++) {
		  for (int i0 = 0; i0 < box_w0; i0++) {
		     buf_ptr[buf_b1+i0] = dat_ptr[dat_b1+i0];
		  }
		  dat_b1 += dat_w0;
		  buf_b1 += box_w0;
	       }
	       dat_b2 += dat_w1 * dat_w0;
	    }
	 }
         break;

      }

      case DATA_TO_VIZ_IS_INT: {

         const tbox::Pointer< pdat::CellData<DIM,int> > ipdata = pdata;
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!ipdata.isNull());
#endif
         const int *const dat_ptr = ipdata->getPointer(depth);

	 if (DIM == 1) {
            int dat_b1 = dat_b2;
	    for (int i0 = 0; i0 < box_w0; i0++) {
	       buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
	    }
	 } else if (DIM == 2) {
            int dat_b1 = dat_b2;
            for (int i1 = 0; i1 < box_w1; i1++) {
               for (int i0 = 0; i0 < box_w0; i0++) {
                  buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
               }
               dat_b1 += dat_w0;
               buf_b1 += box_w0;
            }
	 } else if (DIM == 2) {
	    for (int i2 = 0; i2 < box_w2; i2++) {
	       int dat_b1 = dat_b2;
	       for (int i1 = 0; i1 < box_w1; i1++) {
		  for (int i0 = 0; i0 < box_w0; i0++) {
		     buf_ptr[buf_b1+i0] = (double) dat_ptr[dat_b1+i0];
		  }
		  dat_b1 += dat_w0;
		  buf_b1 += box_w0;
	       }
	       dat_b2 += dat_w1 * dat_w0;
	    }
	 }

         break;

      }

      default: {

      }

   }

}

}
}
#endif
