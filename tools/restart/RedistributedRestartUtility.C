#include "RedistributedRestartUtility.h"

#include "tbox/List.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#include <assert.h>

#define NAME_BUF_SIZE (32)

/*
**************************************************************************
* writeRedistributedRestartFiles                                         *
**************************************************************************
*/

void RedistributedRestartUtility::writeRedistributedRestartFiles(
   const string& output_dirname,
   const string& input_dirname,
   const int total_input_files,
   const int total_output_files,
   const tbox::Array< tbox::Array<int> >& file_mapping,
   const int restore_num)
{
   int num_files_written = 0;
   int num_iterations = tbox::MathUtilities<int>::Min(total_input_files,
                                              total_output_files);

   for (int icount = 0; icount < num_iterations; icount++) {

      //We are writing to one file or reading only one file
      int num_files_to_read = (total_input_files < total_output_files) ?
                              1 : file_mapping[icount].size();
      int num_files_to_write = (total_input_files < total_output_files) ?
                               file_mapping[icount].size() : 1;

      char restore_buf[NAME_BUF_SIZE];
      char nodes_buf[NAME_BUF_SIZE];
      sprintf(restore_buf, "/restore.%05d", restore_num);
      sprintf(nodes_buf, "/nodes.%05d", total_output_files);

      string restart_dirname = output_dirname + restore_buf + nodes_buf;

      //Make the subdirectories if this is the first iteration.
      if (icount == 0) {
         tbox::Utilities::recursiveMkdir(restart_dirname);
      }

      //Mount the output files on an array of output databases
      tbox::Array< tbox::Pointer<tbox::HDFDatabase> >
         output_dbs(num_files_to_write);

      char proc_buf[NAME_BUF_SIZE];
      for (int j = 0; j < num_files_to_write; j++) {

         sprintf(proc_buf, "/proc.%05d", num_files_written+j);

         string output_filename = restart_dirname + proc_buf;

         output_dbs[j] = new tbox::HDFDatabase(output_filename);

         int open_success = output_dbs[j]->mount(output_filename, "W");

         if (open_success < 0) {
            TBOX_ERROR("Failed to open output file " << output_filename <<
                       "  HDF return code:  " << open_success);
         }

      }

      //Mount the input files on an array of input databases.
      tbox::Array< tbox::Pointer<tbox::HDFDatabase> >
         input_dbs(num_files_to_read);

      sprintf(nodes_buf, "/nodes.%05d", total_input_files);

      tbox::Array<string> input_keys(0);
      tbox::Array<string> test_keys(0);
      int num_keys = 0;

      for (int i = 0; i < num_files_to_read; i++) {

         int cur_in_file_id;
         if (total_input_files < total_output_files) {
            //cur_in_file_id = num_files_written;
            cur_in_file_id = icount;
         } else {
            cur_in_file_id = file_mapping[icount][i];
         }  
         sprintf(proc_buf, "/proc.%05d", cur_in_file_id);

         string restart_filename = input_dirname + restore_buf + nodes_buf +
                                   proc_buf;

         input_dbs[i] = new tbox::HDFDatabase(restart_filename);

         int open_success = input_dbs[i]->mount(restart_filename, "R");

         if (open_success < 0) {
            TBOX_ERROR("Failed to open input file " << restart_filename <<
                       "  HDF return code:  " << open_success);
         }

         //Get the array of input keys.
         if (i == 0) {
            input_keys = input_dbs[i]->getAllKeys();
            num_keys = input_keys.size(); 
         } else {
            test_keys = input_dbs[i]->getAllKeys();
            if (test_keys.size() != num_keys) {
               TBOX_ERROR("Input files contain differing number of keys"); 
            }
         }
      }

      //For every input key, call the recursive function that reads from the
      //input databases and writes to output databases.
      for (int j = 0; j < num_keys; j++) {
         readAndWriteRestartData(output_dbs,
                                 input_dbs,
                                 input_keys[j],
                                 &file_mapping,
                                 num_files_written,
                                 icount,
                                 total_input_files,
                                 total_output_files);
      }

      //Unmount the databases.  This closes the files.
      int k;
      for (k = 0; k < num_files_to_read; k++) {
         input_dbs[k]->unmount();
      }
      for (k = 0; k < num_files_to_write; k++) {
         output_dbs[k]->unmount();
      }

      num_files_written += num_files_to_write;
   }
}

/*
**************************************************************************
* readAndWriteRestartData                                                *
**************************************************************************
*/

void RedistributedRestartUtility::readAndWriteRestartData(
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& output_dbs,
   const tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& input_dbs,
   const string& key,
   const tbox::Array< tbox::Array<int> >* file_mapping, // = NULL
   const int num_files_written, // = -1,
   const int input_proc_num, // = -1
   const int total_input_files, // = -1,
   const int total_output_files) // = -1););
{
#ifdef DEBUG_CHECK_ASSERTIONS
   //One of the database arrays must be of size 1, and the other must be of
   //size >= 1.
   assert(output_dbs.size() >= 1);
   assert(input_dbs.size() >= 1);
   assert(input_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //This function works under the assumption that all of the input databases
   //contain the same keys, so we only need to check the type associated with
   //the key with one input database.

   //If the key is associated with any type other than Database, then the data
   //can be read from input_dbs[0] and written to every element of output_dbs.
   //The Database case is handled separately.

   if (input_dbs[0]->isDatabase(key)) {
      tbox::Pointer<tbox::HDFDatabase> db = input_dbs[0]->getDatabase(key);

      if (db->keyExists("d_is_patch_level") &&
          db->getBool("d_is_patch_level")) {

         //Here we are handling the input database(s) for a PatchLevel.

         tbox::Array<int> proc_map(0);

         //A new processor mapping must be created that maps the patches on
         //the level to output files, where each output file represents a
         //processor that will use the redistributed restart output.
         createNewProcessorMapping(proc_map, db, *file_mapping,
                                   total_input_files,
                                   total_output_files);

         //Create array of level input databases.
         tbox::Array< tbox::Pointer<tbox::HDFDatabase> > level_in_dbs;
         level_in_dbs.resizeArray(input_dbs.size());

         for (int i = 0; i < input_dbs.size(); i++) {
            level_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         //input_proc_nums is an array that contains all of the processor
         //numbers that created the input databases that are currently
         //being processed.
         tbox::Array<int> input_proc_nums;
         if (total_input_files < total_output_files) {
            input_proc_nums.resizeArray(1);
            input_proc_nums[0] = input_proc_num;
         } else {
            input_proc_nums = (*file_mapping)[num_files_written];
         }

         //Call routine to write output according to the new processor mapping
         readAndWritePatchLevelRestartData(output_dbs, level_in_dbs,
                                           key, proc_map,
                                           num_files_written,
                                           input_proc_nums);

      } else {

         //If this block is entered, then the key represents a database that
         //is not a patch level database.  We created child database arrays
         //for input_dbs and output_dbs, and then call readAndWriteRestartData
         //recursively.

         tbox::Array< tbox::Pointer<tbox::HDFDatabase> > child_in_dbs;
         child_in_dbs.resizeArray(input_dbs.size());

         for (int i = 0; i < input_dbs.size(); i++) {
            child_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         tbox::Array< tbox::Pointer<tbox::HDFDatabase> > child_out_dbs;
         child_out_dbs.resizeArray(output_dbs.size());

         for (int i = 0; i < output_dbs.size(); i++) {
            child_out_dbs[i] = output_dbs[i]->putDatabase(key);
         }

         tbox::Array<string> child_keys = db->getAllKeys();

         for (int j = 0; j < child_keys.size(); j++) {
            readAndWriteRestartData(child_out_dbs,
                                    child_in_dbs,
                                    child_keys[j],
                                    file_mapping,
                                    num_files_written,
                                    input_proc_num,
                                    total_input_files,
                                    total_output_files);
         }
      }
   } else if (input_dbs[0]->isInteger(key)) {

      tbox::Array<int> int_array(0);

      int_array = input_dbs[0]->getIntegerArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putIntegerArray(key, int_array);
      }

   } else if (input_dbs[0]->isDouble(key)) {

      tbox::Array<double> double_array(0);

      double_array = input_dbs[0]->getDoubleArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putDoubleArray(key, double_array);
      }
 
   } else if (input_dbs[0]->isBool(key)) {

      tbox::Array<bool> bool_array(0);

      bool_array = input_dbs[0]->getBoolArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putBoolArray(key, bool_array);
      }

   } else if (input_dbs[0]->isDatabaseBox(key)) {

      tbox::Array<tbox::DatabaseBox> box_array(0);

      box_array = input_dbs[0]->getDatabaseBoxArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putDatabaseBoxArray(key, box_array);
      }

   } else if (input_dbs[0]->isString(key)) {

      tbox::Array<string> string_array(0);

      string_array = input_dbs[0]->getStringArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putStringArray(key, string_array);
      }

   } else if (input_dbs[0]->isComplex(key)) {

      tbox::Array<dcomplex> complex_array(0);

      complex_array = input_dbs[0]->getComplexArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putComplexArray(key, complex_array);
      }

   } else if (input_dbs[0]->isChar(key)) {

      tbox::Array<char> char_array(0);

      char_array = input_dbs[0]->getCharArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putCharArray(key, char_array);
      }

   } else if (input_dbs[0]->isFloat(key)) {

      tbox::Array<float> float_array(0);

      float_array = input_dbs[0]->getFloatArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putFloatArray(key, float_array);
      }

   } else {

      TBOX_ERROR("The key " << key << " is invalid or not associated with a supported datatype.");

   }
}

/*
**************************************************************************
* readAndWritePatchLevelRestartData                                      *
**************************************************************************
*/

void RedistributedRestartUtility::readAndWritePatchLevelRestartData(
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& output_dbs,
   const tbox::Array< tbox::Pointer<tbox::HDFDatabase> >& level_in_dbs,
   const string& key,
   const tbox::Array<int>& proc_mapping,
   const int num_files_written,
   const tbox::Array<int>& input_proc_nums)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(output_dbs.size() >= 1);
   assert(level_in_dbs.size() >= 1);
   assert(level_in_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //Create an array of level output databases
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> > level_out_dbs;
   level_out_dbs.resizeArray(output_dbs.size());

   for (int i = 0; i < output_dbs.size(); i++) {
      level_out_dbs[i] = output_dbs[i]->putDatabase(key);
   }

   //Read in data that is global to every processor
   bool is_patch_level = level_in_dbs[0]->getBool("d_is_patch_level");
   int version = level_in_dbs[0]->getInteger("HIER_PATCH_LEVEL_VERSION");
   tbox::Array<tbox::DatabaseBox> box_array =
      level_in_dbs[0]->getDatabaseBoxArray("d_boxes");
   int num_patches = level_in_dbs[0]->getInteger("d_number_patches");
   tbox::Array<int> ratio_to_zero =
      level_in_dbs[0]->getIntegerArray("d_ratio_to_level_zero");
   tbox::Array<tbox::DatabaseBox> physical_domain =
      level_in_dbs[0]->getDatabaseBoxArray("d_physical_domain");
   int level_number = level_in_dbs[0]->getInteger("d_level_number");
   int next_coarser_level =
      level_in_dbs[0]->getInteger("d_next_coarser_level_number"); 
   bool in_hierarchy = level_in_dbs[0]->getBool("d_in_hierarchy");
   tbox::Array<int> ratio_to_coarser =
      level_in_dbs[0]->getIntegerArray("d_ratio_to_coarser_level");
   tbox::Array<int> input_mapping =
      level_in_dbs[0]->getIntegerArray("d_mapping");

   //Write out global data.
   for (int i = 0; i < level_out_dbs.size(); i++) {
      level_out_dbs[i]->putBool("d_is_patch_level", is_patch_level);
      level_out_dbs[i]->putInteger("HIER_PATCH_LEVEL_VERSION", version);
      level_out_dbs[i]->putDatabaseBoxArray("d_boxes", box_array);
      level_out_dbs[i]->putInteger("d_number_patches", num_patches);
      level_out_dbs[i]->putIntegerArray("d_ratio_to_level_zero",
                                        ratio_to_zero);
      level_out_dbs[i]->putDatabaseBoxArray("d_physical_domain",
                                           physical_domain);
      level_out_dbs[i]->putInteger("d_level_number", level_number);
      level_out_dbs[i]->putInteger("d_next_coarser_level_number",
                                  next_coarser_level);
      level_out_dbs[i]->putBool("d_in_hierarchy", in_hierarchy);
      level_out_dbs[i]->putIntegerArray("d_ratio_to_coarser_level",
                                       ratio_to_coarser);

      level_out_dbs[i]->putIntegerArray("d_mapping", proc_mapping);
   }

   //Each iteration of this loop processes the patches from one input
   //database.
   for (int i = 0; i < input_proc_nums.size(); i++) {
      int input_id = input_proc_nums[i];

      //This list will contain all of the patch numbers that came from a
      //single processor.
      tbox::List<int> local_patch_nums;

      for (int j = 0; j < num_patches; j++) {
         if (input_mapping[j] == input_id) {
            local_patch_nums.addItem(j);
         }
      }

      //For every patch number, get the patch database from input,
      //create a database for output, and call routine to read and
      //write patch database data.
      for (tbox::List<int>::Iterator it(local_patch_nums); it; it++) {
         int output_id = proc_mapping[it()] - num_files_written;

         char patch_name[NAME_BUF_SIZE];
         sprintf(patch_name, "level_%04d-patch_%04d", level_number, it());

         tbox::Pointer<tbox::HDFDatabase> patch_in_db =
            level_in_dbs[i]->getDatabase(patch_name);

         tbox::Pointer<tbox::HDFDatabase> patch_out_db =
            level_out_dbs[output_id]->putDatabase(patch_name);

         readAndWritePatchRestartData(patch_out_db, patch_in_db);

      }
   }

}

/*
**************************************************************************
* readAndWritePatchRestartData                                           *
**************************************************************************
*/

void RedistributedRestartUtility::readAndWritePatchRestartData(
   tbox::Pointer<tbox::HDFDatabase>& patch_out_db,
   const tbox::Pointer<tbox::HDFDatabase>& patch_in_db)
{
   //Get the keys in the patch input database.
   tbox::Array<string> keys = patch_in_db->getAllKeys();

   //Place the database on arrays of length 1.
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> > in_db_array(1);
   tbox::Array< tbox::Pointer<tbox::HDFDatabase> > out_db_array(1);

   in_db_array[0] = patch_in_db;
   out_db_array[0] = patch_out_db;

   //Call recursive function to read and write the data associated with each
   //key.
   for (int i = 0; i < keys.size(); i++) {
      readAndWriteRestartData(out_db_array, in_db_array, keys[i]);
   }
}

/*
**************************************************************************
* createNewProcessorMapping                                              *
**************************************************************************
*/

void RedistributedRestartUtility::createNewProcessorMapping(
   tbox::Array<int>& proc_map,
   const tbox::Pointer<tbox::HDFDatabase>& level_db,
   const tbox::Array< tbox::Array<int> >& file_mapping,
   const int total_input_files,
   const int total_output_files)
{

   tbox::Array<int> input_proc_map = level_db->getIntegerArray("d_mapping");

   int num_patches = level_db->getInteger("d_number_patches");

   proc_map.resizeArray(num_patches);

   if (total_input_files < total_output_files) {
      for (int i = 0; i < total_input_files; i++) {
         tbox::List<int> input_patch_nums;

         for (int j = 0; j < num_patches; j++) {

            if (input_proc_map[j] == i) {
               input_patch_nums.addItem(j);
            }
         }

         int counter = 0;
         for (tbox::List<int>::Iterator ip(input_patch_nums); ip; ip++) {
            proc_map[ip()] = file_mapping[i][counter];
            counter = (counter+1) % file_mapping[i].size();
         }
      }
   } else {
      for (int i = 0; i < total_input_files; i++) {
         tbox::List<int> input_patch_nums;

         int output_id = 0;
         for (int j = 0; j < total_output_files; j++) {
            int num_inputs = file_mapping[j].size();
            if (file_mapping[j][num_inputs-1] >= i &&
                file_mapping[j][0] <= i) {
               output_id = j;
               break;
            }
         }

         for (int k = 0; k < num_patches; k++) {
            if (input_proc_map[k] == i) {
               input_patch_nums.addItem(k);
            }
         }

         for (tbox::List<int>::Iterator ip(input_patch_nums); ip; ip++) {
            proc_map[ip()] = output_id;
         }
      }
   }
}

