//
// File:        Utilities.C
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A collection of trivial utility functions such as min and max
//

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include "tbox/MPI.h"
#include "tbox/PIO.h"

#include "tbox/Utilities.h"
#ifdef DEBUG_NO_INLINE
#include "tbox/Utilities.I"
#endif

namespace SAMRAI {
   namespace tbox {


bool Utilities::feq(const float a, const float b)
{
   float absmax = fmax(fabs(a),fabs(b));
   return (fabs(a-b)/fmax(absmax, FLT_EPSILON) < sqrt(FLT_EPSILON));
}

bool Utilities::deq(const double a, const double b)
{
   double absmax = dmax(dabs(a),dabs(b));
   return (dabs(a-b)/dmax(absmax, DBL_EPSILON) < sqrt(DBL_EPSILON));
}

void Utilities::renameFile(const char *old_filename,
                                const char *new_filename)
{
   rename(old_filename,new_filename);
}

bool Utilities::ceq(const dcomplex a, const dcomplex b)
{  
   double a_re = real(a);
   double a_im = imag(a);
   double b_re = real(b);
   double b_im = imag(b);

   return ( deq(a_re,b_re) && deq(a_im,b_im) );
}

void Utilities::recursiveMkdir(
   const string& path, 
   mode_t mode,
   bool only_node_zero_creates)
{

#ifdef _MSC_VER
   const char seperator = '/';
#define mkdir(path, mode) mkdir(path)
#else
   const char seperator = '/';
#endif


   if ( (!only_node_zero_creates) || (MPI::getRank() == 0)) {
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
               << "is NOT a directory" << endl);
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
		       << path_buf << endl);
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
			  << path_buf << endl);
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
      MPI::barrier();
   }
}


void Utilities::abort(const string &message, 
		           const string &filename, 
			   const int line) 
{
   perr << "Program abort called in file ``" << filename
        << "'' at line " << line << endl;
   perr << "ERROR MESSAGE: " << endl << message.c_str() << endl;
   perr << flush;

   MPI::abort();
}

void Utilities::printWarning(const string &message, 
        			  const string &filename, 
        			  const int line) 
{
   plog << "Warning in file ``" << filename
        << "'' at line " << line << endl;
   plog << "WARNING MESSAGE: " << endl << message.c_str() << endl;
   plog << flush;
}


}
}

