//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/restart/main.C $
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1806 $
// Modified:    $LastChangedDate: 2007-12-18 22:50:36 -0800 (Tue, 18 Dec 2007) $
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
#include "tbox/SAMRAI_MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAIManager.h"
// Headers for major algorithm/data structure objects


// Header for application-specific algorithm/data structure object


#ifndef NAME_BUFSIZE
#define NAME_BUFSIZE (32)
#endif

using namespace SAMRAI;

#include <iostream>
#include <cstring>
using namespace std;


int main( int argc, char *argv[])
{

   const string slash = "/";

   tbox::SAMRAI_MPI::init(&argc, &argv);
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

   string restore_dirname;
   if( read_dirname.compare(0, 1, "/") == 0 ) {
      restore_dirname = read_dirname + restore_buf;
   } else {
      restore_dirname = "."  + slash + read_dirname + restore_buf;
   }

   if( write_dirname.compare(0, 1, "/") != 0 ) {
      write_dirname = "." + slash + write_dirname;
   }

   struct dirent **namelist;

   //directory should have three entries:  ., .., and nodes.*****
   int num_entries = scandir(restore_dirname.c_str(), &namelist, 0, 0); 

   if(num_entries < 0) {
      TBOX_ERROR("restore directory not found.");
   }

   // Expect only a single run to be in restore directory
   if (num_entries > 3) {
      TBOX_ERROR("restore directory should contain restart files for a single run; this probably indicates runs with different number of nodes have been done in in the restart directory");
   }

   string nodes_dirname;
   string prefix = "nodes.";
   for(int i = 0; i < num_entries; i++) {
      if ( strncmp(prefix.c_str(), namelist[i] -> d_name, prefix.length()) == 0 ) {
	 nodes_dirname = namelist[i]->d_name; 
      }
   }

   int num_input_files = 0;

   // Check if nodes_dirname is valid and extract number of processors for the saved run.
   if ( nodes_dirname.size() == 11 ) {

      string int_str = &(nodes_dirname.c_str()[6]);

      num_input_files = atoi(int_str.c_str());

   } else {
      TBOX_ERROR("nodes.***** subdirectory not found in restore directory.  A directory with a name such as nodes.00016 must be present, with the number indicating the number of process for the run");
   }


   free(namelist);

   string full_nodes_dirname = restore_dirname + slash + nodes_dirname;
   num_entries = scandir(full_nodes_dirname.c_str(), &namelist, 0, 0);
   if (num_entries != num_input_files+2) {
      TBOX_ERROR("number of files in nodes subdirectory does not match the number indicated in the directory's name"); 
   }

   free(namelist);

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
