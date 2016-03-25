//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/Utilities.h $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1768 $
// Modified:    $LastChangedDate: 2007-12-11 16:02:04 -0800 (Tue, 11 Dec 2007) $
// Description: Utility functions for error reporting, file manipulation, etc.
//

#ifndef included_tbox_Utilities
#define included_tbox_Utilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_String
#include <string>
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

/*!
 * Utilities is a Singleton class containing basic routines for error 
 * reporting, file manipulations, etc.
 */

struct Utilities
{
   /*!
    * Create the directory specified by the path string.  Permissions are set 
    * by default to rwx by user.  The intermediate directories in the 
    * path are created if they do not already exist.  When 
    * only_node_zero_creates is true, only node zero creates the 
    * directories.  Otherwise, all nodes create the directories.
    */
   static void recursiveMkdir(const std::string& path, 
			      mode_t mode = (S_IRUSR|S_IWUSR|S_IXUSR),
			      bool only_node_zero_creates = true);

   /*!
    * Rename a file from old file name to new file name.
    */
   static void renameFile(const std::string& old_filename, 
                          const std::string& new_filename);

   /*!
    * Convert an integer to a string.
    *
    * The returned string is padded with zeros as needed so that it
    * contains at least the number of characters indicated by the 
    * minimum width argument.  When the number is positive, the 
    * string is padded on the left. When the number is negative,
    * the '-' sign appears first, followed by the integer value 
    * padded on the left with zeros.  For example, the statement
    * intToString(12, 5) returns "00012" and the statement 
    * intToString(-12, 5) returns "-0012".
    */
   static std::string intToString(int num, int min_width = 1);

   /*!
    * Aborts the run after printing an error message with file and
    * linenumber information.
    */
   static void abort(const std::string &message, 
		     const std::string &filename,
		     const int line);

   /*!
    * Logs warning message with file & location.
    */
   static void printWarning(const std::string &message, 
		            const std::string &filename,
		            const int line);

};


/*!
 * A statement that does nothing, for insure++ make it something 
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if(0) int nullstatement=0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a variable, use to avoid GNU compiler 
 * warnings about unused variables.
 */
#define NULL_USE(variable) do { \
       if(0) {char *temp = (char *)&variable; temp++;} \
    } while (0)

   /*!
    * Throw an error exception from within any C++ source code.  The 
    * macro argument may be any standard ostream expression.  The file and
    * line number of the abort are also printed.
    */
#ifndef LACKS_SSTREAM
#define TBOX_ERROR(X) do {					\
      std::ostringstream tboxos;					\
      tboxos << X << std::ends;					\
      SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_ERROR(X) do {					\
      std::ostrstream tboxos;					\
      tboxos << X << std::ends;					\
      SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /*!
    * Print a warning without exit.  Print file and line number of the warning.
    */
#ifndef LACKS_SSTREAM
#define TBOX_WARNING(X) do {					\
      std::ostringstream tboxos;					\
      tboxos << X << std::ends;					\
      SAMRAI::tbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_WARNING(X) do {					\
      std::ostrstream tboxos;					\
      tboxos << X << std::ends;					\
      SAMRAI::tbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /*!
    * Throw an error exception from within any C++ source code if the
    * given expression is not true.  This is a parallel-friendly version
    * of assert.
    * The file and line number of the abort are also printed.
    */
#ifdef HAVE_STRINGIZE
#ifndef LACKS_SSTREAM
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         std::ostringstream tboxos;                             \
         tboxos << "Failed assertion: " << #EXP << std::ends;        \
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#else
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         std::ostrstream tboxos;                                \
         tboxos << "Failed assertion: " << #EXP << std::ends;        \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#endif
#else
#ifndef LACKS_SSTREAM
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         std::ostringstream tboxos;                             \
         tboxos << "Failed assertion: " << std::ends;                \
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#else
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         std::ostrstream tboxos;                                \
         tboxos << "Failed assertion: " << std::ends;                \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#endif
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
#define PETSC_SAMRAI_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         std::ostringstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#else
#define PETSC_SAMRAI_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         std::ostrstream tboxos;							\
         CHKERRCONTINUE(ierr); 							\
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
      } 									\
} while (0)
#endif
#endif

}
}

#endif
