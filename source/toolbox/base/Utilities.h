//
// File:        Utilities.h
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A collection of trivial utility functions such as min and max
//

#ifndef included_tbox_Utilities
#define included_tbox_Utilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif

#ifndef included_sys_types
#include <sys/types.h>
#define included_sys_types
#endif

#ifndef included_sys_stat
#include <sys/stat.h>
#define included_sys_stat
#endif

namespace SAMRAI {
   namespace tbox {

#ifdef _MSC_VER
#include <sys/types.h>
#include <sys/stat.h>
#include <direct.h>
typedef int mode_t;
#define  S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#define S_IRUSR 0
#define S_IWUSR 0
#define S_IXUSR 0
#endif

/**
 * Class Utilities is a utility that defines simple utility functions 
 * such as min and max.  This class provides a single location of definition
 * for these functions instead of sprinkling them throughout the code.
 */

struct Utilities
{
   /**
    * Calculate the minimum of two integers.
    */
   static int imin(const int a, const int b);

   /**
    * Calculate the maximum of two integers.
    */
   static int imax(const int a, const int b);

   /**
    * Calculate the absolute value of an integer.
    */
   static int iabs(const int a);

   /**
    * Calculate the minimum of two floats.
    */
   static float fmin(const float a, const float b);

   /**
    * Calculate the maximum of two floats.
    */
   static float fmax(const float a, const float b);

   /**
    * Calculate the absolute value of a float.
    */
   static float fabs(const float a);

   /**
    * Returns true if a and b have a relative difference less than
    * sqrt(mach_eps).
    */
   static bool feq(const float a, const float b);

   /**
    * Calculate the minimum of two doubles.
    */
   static double dmin(const double a, const double b);

   /**
    * Calculate the maximum of two doubles.
    */
   static double dmax(const double a, const double b);

   /**
    * Calculate the absolute value of a double.
    */
   static double dabs(const double a);

   /**
    * Returns true if a and b have a relative difference less than
    * sqrt(mach_eps).
    */
   static bool deq(const double a, const double b);

   /**
    * Returns true if both the real and imaginary parts of a and b have a 
    * relative difference less than sqrt(mach_eps).
    */
   static bool ceq(const dcomplex a, const dcomplex b);

   /**
    * Creates the directory specified in path.  Permissions are set 
    * by default to rwx by user.  The intermediate directories in the 
    * path are created if they do not already exist.  When 
    * only_node_zero_creates is true, only node zero creates the 
    * directories.  Otherwise, all nodes create the directories.
    */
   static void recursiveMkdir(const string& path, 
			      mode_t mode = (S_IRUSR|S_IWUSR|S_IXUSR),
			      bool only_node_zero_creates = true);

   /**
    * Rename a file.
    */
   static void renameFile(const char *old_filename, const char *new_filename);

   /**
    * Aborts the run after printing an error message with file and
    * linenumber information.
    */
   static void abort(const string &message, 
		     const string &filename,
		     const int line);

   /**
    * Logs warning message with file & location.
    */
   static void printWarning(const string &message, 
		            const string &filename,
		            const int line);

};


/**
 * A statement that does nothing, for insure++ make it something 
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if(0) int nullstatement=0
#else
#define NULL_STATEMENT
#endif

/**
 * A null use of a variable, use to avoid GNU compiler 
 * warnings about unused variables.
 */
#define NULL_USE(variable) do { \
       if(0) {char *temp = (char *)&variable; temp++;} \
    } while (0)

   /**
    * Throw an error exception from within any C++ source code.  The 
    * macro argument may be any standard ostream expression.  The file and
    * line number of the abort are also printed.
    */
#ifndef LACKS_SSTREAM
#define TBOX_ERROR(X) do {					\
      ostringstream tboxos;					\
      tboxos << X << ends;					\
      SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_ERROR(X) do {					\
      ostrstream tboxos;					\
      tboxos << X << ends;					\
      SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /**
    * Print a warning without exit.  Print file and line number of the warning.
    */
#ifndef LACKS_SSTREAM
#define TBOX_WARNING(X) do {					\
      ostringstream tboxos;					\
      tboxos << X << ends;					\
      SAMRAI::tbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_WARNING(X) do {					\
      ostrstream tboxos;					\
      tboxos << X << ends;					\
      SAMRAI::tbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /**
    * Throw an error exception from within any C++ source code if the
    * given expression is not true.  This is a parallel-friendly version
    * of assert.
    * The file and line number of the abort are also printed.
    */
#ifndef LACKS_SSTREAM
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         ostringstream tboxos;                                  \
         tboxos << "Failed assertion: " << ends;\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#else
#define TBOX_ASSERT(EXP) do {                                    \
      if ( !(EXP) ) {                                            \
         ostrstream tboxos;                                      \
         tboxos << "Failed assertion: " << ends;                 \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
      }                                                          \
} while (0)
#endif


   /**
    * Throw an error exception from within any C++ source code.  This is
    * is similar to TBOX_ERROR(), but is designed to be invoked after a
    * call to a PETSc library function.  In other words, it acts similarly
    * to the PETSc CHKERRQ(ierr) macro.
    */
#ifdef HAVE_PETSC

/*
 * In the following, "CHKERRCONTINUE(ierr);" will cause PETSc to print out
 * a stack trace that led to the error; this may be useful for debugging.
 */

#ifndef LACKS_SSTREAM
/**
 * If using Insure++ then avoid using petsc error macro since it is 
 * a null statement.  Use one that will not cause a warning.
 */
#ifdef __INSURE__
#define PETSC_SAMRAI_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         ostringstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#else
#define PETSC_SAMRAI_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         ostringstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#endif

#else
#define PETSC_SAMRAI_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         ostrstream tboxos;							\
         CHKERRCONTINUE(ierr); 							\
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
      } 									\
} while (0)
#endif
#endif


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Utilities.I"
#endif
#endif
