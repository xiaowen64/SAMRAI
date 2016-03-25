//
// File:	MainRestartData.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Concrete subclass of tbox::Serializable for storing data in main.
//

#ifndef included_MainRestartData
#define included_MainRestartData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

using namespace SAMRAI;

/**
 * Class MainRestartData is a concrete subclass of tbox::Serializable that is 
 * used for storing and accessing the data in main that is necessary for 
 * restart.
 *
 * The input and restart data for the main program are summarized as follows:
 *
 * \verbatim
 * Input:
 *    Required keyword assignments: max_timesteps
 *    Optional keyword assignments: start_time, end_time, regrid_step, 
 *                                  tag_buffer
 *
 * A sample input file entry might look like:
 *
 *   max_timesteps       = 10
 *   start_time          = 0.
 *   end_time            = 10.
 *   regrid_step         = 2
 *   tag_buffer          = 2
 *
 * Restart:
 *    Data written: d_max_timesteps, d_start_time, d_end_time, d_regrid_step,
 *                  d_tag_buffer, d_loop_time, d_iteration_number
 *    Data read:    d_max_timesteps, d_start_time, d_end_time, d_regrid_step,
 *                  d_tag_buffer, d_loop_time, d_iteration_number
 *
 *    Input overwrites all restart values.
 * \endverbatim
 *
 */

class MainRestartData: public tbox::Serializable
{
public:
   /**
    * The constructor for the serializable base class does nothing interesting.
    */
   MainRestartData(const string& object_name,
                   tbox::Pointer<tbox::Database> input_db);

   /**
    * The virtual destructor for the serializable base class does nothing
    * interesting.
    */
   virtual ~MainRestartData();

   /**
    * Returns d_max_timesteps.
    */
   virtual int getMaxTimesteps();

   /**
    * Returns d_start_time.
    */
   virtual double getStartTime();

   /**
    * Returns d_end_time.
    */
   virtual double getEndTime();

   /**
    * Returns d_regrid_step.
    */
   virtual int getRegridStep();

   /**
    * Returns d_tag_buffer.
    */
   virtual int getTagBuffer();

   /**
    * Returns d_loop_time.
    */
   virtual double getLoopTime();

   /**
    * Sets d_loop_time.
    */
   virtual void setLoopTime(const double loop_time);

   /**
    * Returns d_iteration_number.
    */
   virtual int getIterationNumber();

   /**
    * Sets d_iteration_number.
    */
   virtual void setIterationNumber(const int iter_num);

   /**
    * Writes out d_max_timesteps, d_start_time, d_end_time,
    * d_regrid_step, d_tag_buffer, d_loop_time, d_iteration_number.
    */
   virtual void putToDatabase( tbox::Pointer<tbox::Database> db);

private:
   /**
    * Reads in max_timesteps, start_time, end_time,
    * regrid_step, tag_buffer from the specified input database.
    * Any values from the input database override values found
    * in the restart database.
    */
   virtual void getFromInput( tbox::Pointer<tbox::Database> input_db,
                              bool is_from_restart);

   /**
    * Reads in d_max_timesteps, d_start_time, d_end_time,
    * d_regrid_step, d_tag_buffer, d_loop_time, d_iteration_number
    * from the specified restart database.
    */
   virtual void getFromRestart(); 

   int d_max_timesteps;
   double d_start_time;
   double d_end_time;
   int d_regrid_step;
   int d_tag_buffer;
   double d_loop_time;
   int d_iteration_number;

   string d_object_name;
};

#endif
