#include "SAMRAI_config.h"

#ifdef HAVE_HDF5

#include "tbox/Array.h"
#include "tbox/HDFDatabase.h"

#include <string>
using namespace std;
using namespace SAMRAI;

/*
 * RedistributedRestartUtility is a utility class used to build the
 * restart-redistribute tool.  All methods are static and the class
 * contains no data. 
 */

class RedistributedRestartUtility
{

public:

/*
 * Write redistributed restart files to a new restart directory.
 *
 * @param output_dirname      name of directory containing new restart files
 * @param input_dirname       name of directory containing files to be read
 * @param total_input_files   number of input files being read by the tool
 * @param total_output_files  number of output files being created by the tool
 * @param file_mapping        mapping between input and output files
 * @param restore_num         number identifying the restart dump being
 *                            processed
 */
static void writeRedistributedRestartFiles(
   const string& output_dirname,
   const string& input_dirname,
   const int total_input_files,
   const int total_output_files,
   const tbox::Array< tbox::Array<int> >& file_mapping,
   const int restore_num);

private:

/*
 * A recursive function that does the reading and writing of data.  If the
 * key represents any type other than Database, data is read from input and
 * immediately written to output.  If it is a Database, this function is
 * called recursively on every key in that Database.  If the Database
 * represents a PatchLevel, readAndWritePatchLevelRestartData is called
 * instead, in order to properly distribute patches among the output files.
 *
 * The optional arguments are not used when this function is called for
 * reading and writing data from a patch, because the data for a patch
 * exists in only one input file and is written to only one output file.
 * Meaningful values should be given for all of the optional arguments in
 * all other cases.
 */
static void readAndWriteRestartData(
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& output_dbs,
   const tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& input_dbs,
   const string& key,
   const tbox::Array< tbox::Array<int> >* file_mapping = NULL,
   int num_files_written = -1,
   int input_proc_num = -1,
   int total_input_files = -1,
   int total_output_files = -1);

/*
 * Creates a new processor mapping to map patches on a level to the
 * processors represented by the output files.  New mapping is stored
 * in proc_map.
 */
static void createNewProcessorMapping(
   tbox::Array<int>& proc_map,
   const tbox::Pointer<tbox::HDFDatabase>& level_db,
   const tbox::Array< tbox::Array<int> >& file_mapping,
   const int total_input_files,
   const int total_output_files);

/*
 * Reads and writes data in a Database that represents a PatchLevel.
 * Global data are written to every output database, and data representing
 * patches are written to only one output database, determined by
 * proc_mapping.
 */
static void readAndWritePatchLevelRestartData(
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& output_dbs,
   const tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& level_in_dbs,
   const string& key,
   const tbox::Array<int>& proc_mapping,
   const int num_files_written,
   const tbox::Array<int>& input_proc_nums);

/*
 * Reads patch restart data from input database and writes it to output
 * database.
 */
static void readAndWritePatchRestartData(
   tbox::Pointer<tbox::HDFDatabase>& patch_out_db,
   const tbox::Pointer<tbox::HDFDatabase>& patch_in_db);

};

#endif
