//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/Utilities.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2039 $
// Modified:    $LastChangedDate: 2008-03-11 13:23:52 -0700 (Tue, 11 Mar 2008) $
// Description: Utility functions for error reporting, file manipulation, etc.
//

#include "tbox/Utilities.h"

#include "tbox/SAMRAI_MPI.h"
#include "tbox/Logger.h"
#include "tbox/PIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

namespace SAMRAI {
   namespace tbox {

/*
 * Routine to rename a file.
 */
void Utilities::renameFile(const std::string& old_filename,
                           const std::string& new_filename)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!old_filename.empty());
   TBOX_ASSERT(!new_filename.empty());
#endif
   rename(old_filename.c_str(), new_filename.c_str());
}

/*
 * Routine to recursively construct directories based on a relative path name.
 */
void Utilities::recursiveMkdir(
   const std::string& path, 
   mode_t mode,
   bool only_node_zero_creates)
{

#ifdef _MSC_VER
   const char seperator = '/';
#define mkdir(path, mode) mkdir(path)
#else
   const char seperator = '/';
#endif


   if ( (!only_node_zero_creates) || (SAMRAI_MPI::getRank() == 0)) {
      int length = path.length();
      char *path_buf= new char[length+1];
      sprintf(path_buf,"%s",path.c_str());
      struct stat status;
      int pos = length - 1;
   
      /* find part of path that has not yet been created */
      while ( (stat(path_buf,&status) != 0) && (pos >= 0) ) {
   
         /* slide backwards in string until next slash found */
         bool slash_found = false;
         while ( (!slash_found) && (pos >= 0) ) {
           if (path_buf[pos] == seperator) {
              slash_found = true;
              if (pos >= 0) path_buf[pos] = '\0';
           } else pos--;
         }
      } 

      /* 
       * if there is a part of the path that already exists make sure
       * it is really a directory
       */
      if (pos >= 0) {
         if ( !S_ISDIR(status.st_mode) ) {
            TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
               << "    Cannot create directories in path = " << path
               << "\n    because some intermediate item in path exists and"
               << "is NOT a directory" << std::endl);
         }
      }
   
      /* make all directories that do not already exist */
   
      /* 
       * if (pos < 0), then there is no part of the path that
       * already exists.  Need to make the first part of the 
       * path before sliding along path_buf.
       */
      if (pos < 0) {
	 if(mkdir(path_buf,mode) != 0) {
	    TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
		       << "    Cannot create directory  = " 
		       << path_buf << std::endl);
	 }
	 pos = 0;
      }
   
      /* make rest of directories */
      do {
   
         /* slide forward in string until next '\0' found */
         bool null_found = false;
         while ( (!null_found) && (pos < length) ) {
           if (path_buf[pos] == '\0') {
              null_found = true;
              path_buf[pos] = seperator;
           }
           pos++;
         }
   
         /* make directory if not at end of path */
	 if (pos < length) {
	    if(mkdir(path_buf,mode) != 0) {
	       TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
			  << "    Cannot create directory  = " 
			  << path_buf << std::endl);
	    }
	 }
      } while (pos < length);

      delete [] path_buf;
   }

   /* 
    * Make sure all processors wait until node zero creates 
    * the directory structure.
    */
   if (only_node_zero_creates) {
      SAMRAI_MPI::barrier();
   }
}

/*
 * Routine to convert an integer to a string. 
 */
std::string Utilities::intToString(int num, int min_width)
{
   int tmp_width = ( min_width > 0 ? min_width : 1 );
   std::ostringstream os;
   if ( num < 0 ) {
      os << '-' << std::setw(tmp_width-1) << std::setfill('0') << -num;
   } else {
      os << std::setw(tmp_width) << std::setfill('0') << num;
   }
   os << std::flush;

   return(os.str()); //returns the string form of the stringstream object
}

std::string Utilities::nodeToString(int num) 
{
   return intToString(num, 5);
}

std::string Utilities::processorToString(int num) 
{
   return intToString(num, 5);
}


std::string Utilities::patchToString(int num) {
   return intToString(num, 4);
}


std::string Utilities::levelToString(int num) {
   return intToString(num, 4);
}

std::string Utilities::blockToString(int num) {
   return intToString(num, 4);
}

/*
 * Routine that calls abort and prints calling location to error stream.
 */
void Utilities::abort(const std::string &message, 
	              const std::string &filename, 
	              const int line) 
{
   tbox::Logger::getInstance() -> logAbort(message, filename, line);

   SAMRAI_MPI::abort();
}

}
}

