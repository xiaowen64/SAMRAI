//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/plotting/CartesianVizamraiDataWriter.h $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Simple tool to facilitate dumping data to file for Vizamrai
//

#ifndef included_appu_CartesianVizamraiDataWriter
#define included_appu_CartesianVizamraiDataWriter

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "CartesianPatchGeometry.h"
#include "BoxList.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "VisDerivedDataStrategy.h"
#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "tbox/FileStream.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace appu {

/*!
 * Class CartesianVizamraiDataWriter<DIM> is used by SAMRAI-based 
 * application codes to generate Vizamrai data files.  Vizamrai provides 
 * a wide range of visualization and post-processing capabilities 
 * for data generated on a structured AMR patch hierarchy.  Vizamrai 
 * supports cell-centered mesh data.  This class supports cell-centered
 * data where the underlying data type is either double, float, or int. 
 *
 * This data file writer class can be used to generate Vizamrai data
 * files when the underlying mesh geometry is managed by a
 * geom::CartesianGridGeometry<DIM> object.  Two kinds of data quantities are 
 * supported by this class.  The first is data that resides on an AMR patch 
 * hierarchy when a data file is generated.  The second is "derived" data 
 * that can be computed using data on an AMR patch hierarchy.  The first data 
 * type requires no user intervention when generating a plot file.  The second 
 * type, or "derived" data, requires that a user implement a concrete class 
 * derived from the VisDerivedDataStrategy<DIM> abstract interface.  In 
 * an overloaded function of that class, a user generates the "derived" data 
 * quantity and writes it to a specified file stream.   In either case, scalar or
 * variable data quantities are supported.
 *
 * Before a Vizamrai data writer object can be used to generate Vizamrai
 * data files, it must first be constructed and initialized.  Initialization
 * involves things like setting information about the levels in the AMR 
 * hierarchy that will be plotted, setting the type of the data to write
 * out (e.g., double or float), setting the directory into which plot
 * files are written, and registering variable quantities to send
 * to the plot file.  After initialization, the object can be used to 
 * generate a series of data files during the execution of some simulation
 * code.
 *
 * Typical usage of this Vizamrai data writer object involves several
 * steps among which are variations.  Basic usage is as follows:
 * -# Create a Cartesian Vizamrai data writer object.
 * -# Set the finest level to plot with the member function 
 *       setFinestLevelToPlot().  This is the finest level allowable 
 *       in any plot file; however, only levels in the hierarchy will
 *       be written out.
 * -# For each level OTHER THAN THE COARSEST , specify
 *       the ratio between the index space of the level and the next
 *       coarser level in the AMR hierarchy.  This information 
 *       is needed for scaling the AMR levels properly in Vizamrai.  The
 *       member function to use is setRatioToCoarserLevel().
 * -# Set the type of data to put in the Vizamrai file, using either
 *       setPlotDataToDouble(), or setPlotDataToFloat().  The default is
 *       to write float data.
 * -# Set the directory for the plot files if desired using the function
 *       setDirectoryName().  If no directory is specified, files will
 *       be written into the current directory. 
 * -# If any "derived" quantities will be generated, set the concrete
 *       derived data writer object using setDerivedDataWriter().  
 *       Alternatively, one may specify a derived data writer when registering
 *       each derived data quantity.  This provides flexibility in using 
 *       different derived data writers for different quantities if desired.
 * -# Register data quantities using either registerPlotScalar()/registerPlotVector() or
 *       registerDerivedPlotScalar()/registerDerivedPlotVector().  Derived data quantities 
 *       require only a string identifier and a derived data writer if not specified
 *       earlier.  All other quantities require a string identifier and
 *       an index into the patch data array on the AMR hierarchy.  A scale
 *       factor may also be specified.  Note that you cannot register more
 *       than one variable with the same name.  If this is attempted, an 
 *       error message is logged and the program will exit.
 * -# The resetLevelPlotScalar()/resetLevelPlotVector() member functions 
 *       are provided for cases when a plot quantity lives at different 
 *       patch data indices on different levels.  Calling this routine 
 *       redefines the patch data index for a given quantity on a given 
 *       level that will be written to a plot file.  Before this function 
 *       is called, the quantity must be registered using the 
 *       registerPlotScalar()/registerPlotVector() function.
 * -# Finally, the writer can be used to generate Vizamrai data files
 *       using the member function writePlotData().  Minimally, only
 *       a hierarchy and a file name is needed.  One can also supply an
 *       integer file extension and a plot time, both of which are useful 
 *       when generating visualization files for time-dependent applications.
 *
 * During the course of generating Vizamrai data files, certain features
 * of the plot files may change.  For example, the directory name can be
 * reset using the setDirectoryName() function.  Also, the finest level
 * to plot using the setFinestLevelToPlot() function.  Some restrictions
 * apply.  Please consult the member function documentation for more details.
 * 
 * @see appu::VisDerivedDataStrategy
 */

template<int DIM> class CartesianVizamraiDataWriter : public virtual tbox::DescribedClass
{
public:
   /*!
    * The constructor for CartesianVizamraiDataWriter<DIM> initializes
    * this Vizamrai data writer to a default state.  The default state is
    * to write data as floats (instead of doubles).  Before it can be
    * used for writing data to a Vizamrai plot file, the maximum number
    * of hierarchy levels must be set and the variables to plot must be
    * registered.  See the member functions setMaxHierarchyLevels() and
    * the register functions.  
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty.  This name is used primarily for error reporting.
    *
    * @param name   Name given to object.
    */
   CartesianVizamraiDataWriter(const std::string& name);

   /*!
    * The destructor for the writer does nothing interesting.
    */
   virtual ~CartesianVizamraiDataWriter<DIM>();

   /*!
    * Set the finest hierarchy level to write to a Vizamrai file. 
    * Calling this function is optional.  If this function is not called, 
    * the finest level to be plotted is the finest level in the hierarchy
    * when the writePlotData() function is called or the finest level for 
    * which the function setRatioToCoarsrLevel() has been called, whichever
    * is coarser.  Thus, this routine can be used to increase or decrease 
    * the number of levels to plot, but it cannot be used to plot a level 
    * for which scaling ratio information has not been previously provided.  
    *
    * Calling this function a negative value will result in a plot file 
    * with no data written in it when writePlotData() is called.
    *
    * If the integer is negative, no levels will be written to any plot 
    * file.  If an attempt is made to set this value larger than the finest 
    * level for which previous scaling information has been given (e.g., 
    * using setRatioToCoarserLevel()), a warning message will be logged.
    */
   void setFinestLevelToPlot(int finest_level_number);

   /*!
    * Set the vector ratio between the index space of the patch level with
    * the given number and the next coarser patch level.  This information 
    * is needed for proper scaling of the domain when Vizamrai visualizes 
    * datasets.  
    * 
    * If this Vizamrai data writer object will be used to write multiple
    * data files (e.g., for a series of timesteps) and the number of levels
    * will never change, then calling this routine is optional.  However,
    * if the number of levels to plot may change from one plot file to the 
    * next, this routine must be called for each level finer than the 
    * finest level written in the first Vizamrai plot file written by this 
    * writer object.  
    *
    * In any case, the number of levels actually written to a plot file 
    * will depend on the number of levels in the patch hierarchy sent to 
    * the writePlotData() routine, the finest level indicated in the most
    * recent call to the function setFinestLevelToPlot(), and the manner
    * in which this routine is called.  Also, when used, this function 
    * must be called before the first plot file is written; i.e., before 
    * the writePlotData() member function is called.
    * 
    * When assertion checking is active, a non-recoverable exception results
    * if level number is negative.  If an attempt is made to call this 
    * function after a Vizamrai file has been written by this writer, 
    * an error message will be logged and the program will terminate.
    */
   void setRatioToCoarserLevel(int level_number,
                               const hier::IntVector<DIM>& ratio);

   /*!
    * Set the format of the plot data to type double.   The default 
    * data format is float.  If an attempt is made to change the data
    * format after a Vizamrai file has been written by this writer, a 
    * warning message will be logged.
    */
   void setPlotDataToDouble();

   /*!
    * Set the format of the plot data to type float.   The default
    * data format is float.  If an attempt is made to change the data
    * format after a Vizamrai file has been written by this writer, a
    * warning message will be logged.
    */
   void setPlotDataToFloat();

   /*!
    * Set default user-defined data writer for each derived variable
    * quantity to write to the plot file.  If a non-null derived 
    * data writer is supplied to the member function 
    * registerDerivedPlotScalar()/registerDerivedPlotVecotr(), it will supercede 
    * the one given here.
    * 
    * When assertion checking is active, a non-recoverable exception results
    * if the pointer is null. 
    */
   void setDerivedDataWriter(VisDerivedDataStrategy<DIM>* derived_writer);

   /*!
    * Set name of directory into which Vizamrai files will be written.
    * Permissions are set by default to rwx by user.  Any intermediate 
    * directories (including the top level) in the path are created 
    * if they do not already exist.  When running jobs in parallel, only
    * node zero will create the directories.  Calling this function is
    * optional.  If no directory name is specified, the Vizamrai file 
    * will be written into the current directory (i.e., the directory 
    * in which the code is executed).
    * 
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty.
    */
   void setDirectoryName(const std::string& directory_name);

   /*!
    * Register a non-derived scalar variable quantity with given name and patch 
    * data id with this Vizamrai data writer.  The integer data identifier indicates 
    * the descriptor index at which the data may be found on each patch in the
    * hierarchy.  The depth index refers to the depth into the patch data 
    * array of the scalar plot quantity.  For example, when plotting individual 
    * entries of a quantity that is maintained as a vector in each cell on
    * the hierarchy, this is the index into the vector (e.g., think 
    * "y"-component of velocity, where the velocity is stored as a patch
    * data type of depth DIM.  If the data id integer does not correspond 
    * to a valid patch data index, or if the depth id integer is out of
    * range for the given data id, an error message will result and the 
    * program will exit.  See the function registerPlotVector() below for plotting
    * a vector quantity as a vector rather than multiple scalar quantities. 
    *
    * If the scale factor is specified, each data value will be multiplied
    * by this factor before it is written to the data file.
    *
    * The integer data identifier will be checked to see if it is a valid
    * cell-centered patch data object in the patch descriptor.  Valid cell-
    * centered types are double, float, and int.
    *
    * If a variable was previously registered with the same name, the most
    * recent registration replaces the previous one.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty, or if the data id or depth id are negative.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param data_id        Cell-centered patch data descriptor identifyer for 
    *                 the variable (may be double, float, or int).
    * @param depth_id       Depth into the patch data array of the scalar plot
    *                 quantity.
    * @param scale_factor   Scaling factor applied to the variable data when writing
    *                 viz files.
    */
   void registerPlotScalar(const std::string& variable_name,
                           int data_id,
                           int depth_id = 0,
                           double scale_factor = 1.0);

   /*!
    * Change a previously-registered, non-derived variable quantity with 
    * given name with this Vizamrai data writer.  The change redefines
    * the patch data objects written to the plot file on the specified level
    * to the data at the given descriptor index and depth index.  This 
    * function is used when a particular data quantity lives at different
    * patch data slots on different hierarchy levels.  For example, suppose 
    * a plot data quantity lives at a patch data index on every level except
    * the finest hierarchy level, where it lives at a different index.  First,
    * that quantity must be registered using registerPlotScalar()/registerPlotVector().  
    * Second, the patch data index for the finest hierarchy level is reset using this
    * function.  When the data is plotted, it will appear on all levels in
    * the hierarchy.
    *
    * Before this function can be called, the plot quantity must already
    * be registered using a registerPlotScalar() function.  
    * If it is not, then a warning message will result and the variable will
    * not be registered for plotting. 
    *
    * The data id integer argument must correspond to a valid patch 
    * data index with the same type as the data for which the quantity was 
    * originally registered using registerPlotScalar().  
    * The depth id integer (if given) only applies to the case where the original
    * registration of the plot quantity with the string name was as a scalar.  In this
    * case, the depth id must be valid for the given data id.  If either 
    * of these arguments is illegal, an error message will be printed
    * and the program will abort.  See comments for the registerPlotScalar() 
    * function for restrictions on the data types allowed.
    * 
    * If a variable was previously registered with the same name, the most
    * recent registration replaces the previous one.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty, or if the data id, depth id, or level number
    * are negative.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param level_number   The level number on which variable data is defined.
    * @param data_id        Cell-centered patch data descriptor identifyer for 
    *                 the variable (may be double, float, or int).
    * @param depth_id       Depth into the patch data array of the scalar plot
    *                 quantity.
    */
   void resetLevelPlotScalar(const std::string& variable_name,
                               int level_number,
                               int data_id,
                               int depth_id = 0);

   /*!
    * Register a non-derived vector quantity with given name and patch data id 
    * with this Vizamrai data writer.  The integer data identifier indicates the
    * descriptor index at which the data may be found on each patch in the
    * hierarchy.  If the data id integer does not correspond 
    * to a valid patch data index, an error message will result and the 
    * program will exit.  No depth index is required, since data at all depths will 
    * be written for a vector quantitiy.  See the function registerPlotVariable() above 
    * for plotting components of a vector as multiple scalar quantities.
    *
    * If the scale factor is specified, each data value will be multiplied
    * by this factor before it is written to the data file.
    *
    * The integer data identifier will be checked to see if it is a valid
    * cell-centered patch data object in the patch descriptor.  Valid cell-
    * centered types are double, float, and int.
    *
    * If a variable was previously registered with the same name, the most
    * recent registration replaces the previous one.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty, or if the data id is negative.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param data_id        Cell-centered patch data descriptor identifyer for 
    *                 the variable (may be double, float, or int).
    * @param scale_factor   Scaling factor applied to the variable data when writing
    *                 viz files.
    */
   void registerPlotVector(const std::string& variable_name,
                           int data_id,
                           double scale_factor = 1.0);

   /*!
    * Change a previously-registered, non-derived variable quantity with 
    * given name with this Vizamrai data writer.  The change redefines
    * the patch data objects written to the plot file on the specified level
    * to the data at the given descriptor index and depth index.  This 
    * function is used when a particular data quantity lives at different
    * patch data slots on different hierarchy levels.  For example, suppose 
    * a plot data quantity lives at a patch data index on every level except
    * the finest hierarchy level, where it lives at a different index.  First,
    * that quantity must be registered using registerPlotVector().  
    * Second, the patch data index for the finest hierarchy level is reset using this
    * function.  When the data is plotted, it will appear on all levels in
    * the hierarchy.
    *
    * Before this function can be called, the plot quantity must already
    * be registered using a registerPlotVector() function.  
    * If it is not, then a warning message will result and the variable will
    * not be registered for plotting. 
    *
    * The data id integer argument must correspond to a valid patch 
    * data index with the same type as the data for which the quantity was 
    * originally registered using registerPlotVector().  
    * 
    * If a variable was previously registered with the same name, the most
    * recent registration replaces the previous one.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the string is empty, or if the data id, depth id, or level number
    * are negative.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param level_number   The level number on which variable data is defined.
    * @param data_id        Cell-centered patch data descriptor identifyer for 
    *                       the variable (may be double, float, or int).
    */
   void resetLevelPlotVector(const std::string& variable_name,
                             int level_number,
                             int data_id); 
  
   /*!
    * Register a scalar derived variable quantity with the given name with this
    * Vizamrai data writer.  If a non-null derived data writer is given
    * it will supercede any that has previously been specified in the 
    * setDerivedDataWriter() member function.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the argument string is empty.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param derived_writer tbox::Pointer to the VizamraiDerivedDataStrategy class
    *                 (which implements 'writeDerivedDataToStream').
    */
   void registerDerivedPlotScalar(
      const std::string& variable_name,
      VisDerivedDataStrategy<DIM>* derived_writer = 
         (VisDerivedDataStrategy<DIM>*)NULL);

   /*!
    * Register a vector derived variable quantity with given name and 
    * depth (i.e., number of components) with this Vizamrai data writer.  
    * If a non-null derived data writer is given it will supercede any that 
    * has previously been specified in the setDerivedDataWriter() member function.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the argument string is empty.
    *
    * @param variable_name  Vizamrai name identifier for the variable.
    * @param depth          number of components
    * @param derived_writer tbox::Pointer to the VizamraiDerivedDataStrategy class
    *                 (which implements 'writeDerivedDataToStream').
    */
   void registerDerivedPlotVector(
      const std::string& variable_name,
      int depth,
      VisDerivedDataStrategy<DIM>* derived_writer =
         (VisDerivedDataStrategy<DIM>*)NULL);  

   /*!
    * Write Vizamrai plot file with given file name, using data in given
    * hierarchy that matches registered plot quantity information.  The
    * data file will be written in the directory specified in the most
    * recent call to setDirectoryName().  If no directory name has been 
    * specified, the Vizamrai file will be written into the current 
    * directory (i.e., the directory in which the code is executed).
    * Optional arguments are available to provide an integer file number
    * extension and a plot time value to stamp the data in the Vizamrai 
    * plot file.  If a integer extension (>= 0) is provided, the name of the
    * file will contain `file_name.extention'.  Otherwise, the extension
    * will be ignored.  If the time is not specified, a default value of 
    * zero is used.
    *
    * If a variable was previously registered with the same name, the most
    * recent registration replaces the previous one.
    *
    * When assertion checking is active, a non-recoverable exception results
    * if the file name string is empty, or if the hierarchy pointer is null.
    *
    * @param hierarchy hier::Patch hierarchy on which data is defined.
    * @param file_name File prefix for vizamrai data files 
    * @param extension Extension which may be appended to the file name
    *            (i.e. 'file_name.extension'
    * @param plot_time Simulation time of the plot
    */
   void writePlotData(
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const std::string& file_name,
      int extension = -1,
      double plot_time = 0.0);

   /*!
    * Send all data in this Vizamrai file writer to given output stream.
    */
   virtual void printClassData(std::ostream& os) const;

private:
   // The following two members are not implemented
   CartesianVizamraiDataWriter(
      const CartesianVizamraiDataWriter<DIM>&); 
   void operator=(const CartesianVizamraiDataWriter<DIM>&);

   /*
    * The following structure is used to store data about each item
    * to be written to a plot file.
    */
   template<int DIMENSION> struct VizamraiItem {
      std::string d_variable_name;     // name of variable in viz file
      int d_data_id;              // master descriptor id for plot data
      tbox::Array<int> d_level_data_id;  // desc ids for plot data by level
      bool d_isa_vector;          // flag indicating if this item is to 
                                  // be written as a vector
      int d_depth;                // vector depth
      int d_depth_id;             // master depth index for plot data
      tbox::Array<int> d_level_depth_id; // depth indices for plot data by level
      double d_scale_factor;    // scale factor for plotting
      VisDerivedDataStrategy<DIMENSION>* d_derived_writer;
                                // non-NULL if data is derived; NULL otherwise.
      int    d_data_type;       // type of variable data (not necessarily
                                // the same type that is written to file)
   };


   /*
    * Utility routine to register standard (non-derived) plot quantity
    * (either vector or single value variable).
    */
   void registerPlotItem(
      const std::string& variable_name, 
      int data_id,
      bool isa_vector,
      int depth_id,
      double scale_factor);

   /*
    * Utility routine to register derived plot quantity
    * (either vector or single value variable).
    */
   void registerDerivedPlotItem(
      const std::string& variable_name,
      int depth,
      VisDerivedDataStrategy<DIM>* derived_writer);

   /*
    * Return string representing complete name of output file given file name, directory 
    * name, plot file counter, and processor number.
    */
   std::string getFileStreamName(
      const std::string& file_name,
      int extention,
      bool istemporaryflag = false) const;

   /*
    * Utility function that writes Vizamrai header file information.
    */
   void writeVizamraiHeaderInfoToFile(
      tbox::FileStream& file,
      double plot_time,
      int num_outboxes) const;

   /*
    * Utility function to compute a list of plot box regions on
    * each hierarchy level and return the total number of boxes.
    */
   int computeOutputBoxes(
      tbox::Array< hier::BoxList<DIM> > outboxes[], 
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy, 
      int coarsest_plot_level, 
      int finest_plot_level) const;

   /*
    * Utility function to write boundaries of hierarchy patches to 
    * file stream.
    */
   void writePatchBoundariesToFile(
      tbox::FileStream& file,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int coarsest_plot_level,
      int finest_plot_level) const;

   /*
    * Utility function to write variable data to plot file.
    */
   void writeVizamraiVariablesToFile(
      tbox::FileStream& file,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int coarsest_plot_level,
      int finest_plot_level,
      const tbox::Array< hier::BoxList<DIM> > outboxes[]) const;

   /*
    * Utility function to write single plot box geometry data to file.
    */
   void writePlotBoxDataToFile(
      tbox::FileStream& file,
      int level_number,
      const double* const xlo_domain,
      const double* const dx_level,
      const hier::Box<DIM>& plot_box) const;

   /*
    * Utility function to check that information, such as scaling factors 
    * between levels to plot, has been properly specified.  Return true 
    * if data is good and plot file can be written; otherwise return false.
    */
   bool checkLevelInformation(
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int coarsest_plot_level,
      int finest_plot_level);

   /*
    * Utility functions to pack plot data into buffers.
    */
   void packPatchDataIntoDoubleBuffer(
      tbox::Pointer< hier::PatchData<DIM> > pdata,
      int depth,
      int type_of_data,
      const hier::Box<DIM>& box,
      double* buffer) const;

   /*
    * Name of this Vizamrai data writer object (passed into constructor)
    */
   std::string d_object_name;

   /*
    * Boolean flag used to know whether this object has successfully 
    * written a plot file.
    */
   bool d_wrote_plot_file;

   /*
    * Finest hierarchy level information to write to a file.
    */
   int d_finest_plot_level_in_all;
   int d_finest_plot_level_now;

   /*
    * tbox::Array of mesh-scaling ratios from each level to reference level
    * (i.e., coarsest level).
    */
   tbox::Array< hier::IntVector<DIM> > d_scaling_ratios;

   /*
    * Value used to scale "extra" dimensions of plot domain when 
    * DIM < 3.   It is computed when a Vizamrai file is generated.
    * (i.e., coarsest level).
    */
   double d_domain_scale_length;

   /*
    * Type of data to write to Vizamrai plot file.
    * 0 indicates double data, 1 indicates float data (default).
    */
   int d_plot_type;

   /*
    * Default data writer for user defined data.
    */
   VisDerivedDataStrategy<DIM>* d_default_derived_writer;

   /*
    * Default directory into which Vizamrai files will be written.
    */
   std::string d_directory_name;

   /*
    * Master list of variable items written to plot files.
    */
   tbox::List<VizamraiItem<DIM> > d_plot_items;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CartesianVizamraiDataWriter.C"
#endif
