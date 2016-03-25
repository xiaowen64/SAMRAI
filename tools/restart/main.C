//
// File:        main.C
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Main program restart-redistribute tool.
//

#include "SAMRAI_config.h"

#include "RedistributedRestartUtility.h"
#include "tbox/HDFDatabase.h"
#include "tbox/RestartManager.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <unistd.h>
using namespace std;

// Headers for basic SAMRAI objects

#include "tbox/HDFDatabase.h"
#include "tbox/MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAIManager.h"
// Headers for major algorithm/data structure objects


// Header for application-specific algorithm/data structure object


#ifndef NAME_BUFSIZE
#define NAME_BUFSIZE (32)
#endif

using namespace SAMRAI;

int main( int argc, char *argv[])
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   string read_dirname;
   string write_dirname;
   int restore_num = 0;
   int num_output_files = 1;

   if ( (argc != 5) ) {
	 tbox::pout << "USAGE:  " << argv[0] << " input-dir "
	      << "output-dir restore-number num-output-files\n"
	      << endl;
	 exit(-1);
	 return (-1);
   } else {
      read_dirname = argv[1];
      write_dirname = argv[2];
      restore_num = atoi(argv[3]);
      num_output_files = atoi(argv[4]);
   }

   // make string "input-dir/restore.*****
   char restore_buf[NAME_BUFSIZE];

   sprintf(restore_buf,"/restore.%05d",restore_num);

   char working_dir[NAME_BUFSIZE*10];

   getcwd(working_dir, NAME_BUFSIZE*10);

   string slash = "/";
   string restore_dirname = working_dir + slash + read_dirname + restore_buf;

   struct dirent **namelist;

   //directory should have three entries:  ., .., and nodes.*****
   int num_entries = scandir(restore_dirname.c_str(), &namelist, 0, 0); 

   if (num_entries != 3) {
      TBOX_ERROR("restore directory should contain only nodes subdirectory.");
   }

   string nodes_dirname = namelist[num_entries-1]->d_name; 
   int num_input_files = 0;

   //reading in directory name.
   if ((nodes_dirname.size() == 11) && (nodes_dirname.find('n',0) == 0) &&
       (nodes_dirname.find('o',0) == 1) && (nodes_dirname.find('d',0) == 2) &&
       (nodes_dirname.find('e',0) == 3) && (nodes_dirname.find('s',0) == 4) &&
       (nodes_dirname.find('.',0) == 5)) {

      string int_str = &(nodes_dirname.c_str()[6]);

      num_input_files = atoi(int_str.c_str());

   } else {
      TBOX_ERROR("nodes.***** subdirectory not found in restore directory.  A directory with a name such as nodes.00016 must be present, with the number indicating the number of restart files contained within");
   }

   string full_nodes_dirname = restore_dirname + slash + nodes_dirname;
   num_entries = scandir(restore_dirname.c_str(), &namelist, 0, 0);
   if (num_entries != num_input_files+2) {
      TBOX_ERROR("number of files in nodes subdirectory does not match the number indicated in the directory's name"); 
   }

   // file_mapping will have size equal to the lesser value of
   // num_input_files and num_output_files.
   tbox::Array< tbox::Array<int> > file_mapping;
   int file_ratio;
   int remainder;
   if (num_output_files > num_input_files) {
      file_mapping.resizeArray(num_input_files);
      file_ratio = num_output_files / num_input_files;
      remainder = num_output_files % num_input_files;
   } else {
      file_mapping.resizeArray(num_output_files);
      file_ratio = num_input_files / num_output_files;
      remainder = num_input_files % num_output_files;
   }

   int file_counter = 0;

   // fill file_mapping.
   int i;
   for (i = 0; i < file_mapping.size(); i++) {
      if (i < remainder) {
         file_mapping[i].resizeArray(file_ratio+1);
         for (int j = 0; j <= file_ratio; j++) {
            file_mapping[i][j] = file_counter + j;
         }
         file_counter += file_ratio + 1;
      } else {
         file_mapping[i].resizeArray(file_ratio);
         for (int j = 0; j < file_ratio; j++) {
            file_mapping[i][j] = file_counter + j;
         }
         file_counter += file_ratio;
      }
   }

   //Write all of the redistributed restart files.
   RedistributedRestartUtility::writeRedistributedRestartFiles(
      write_dirname,
      read_dirname,
      num_input_files,
      num_output_files,
      file_mapping,
      restore_num);

   tbox::SAMRAIManager::shutdown();

   return(0);

}
