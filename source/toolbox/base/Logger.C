//
// File:        $URL$
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2028 $
// Modified:    $LastChangedDate: 2008-02-29 13:26:00 -0800 (Fri, 29 Feb 2008) $
// Description: Utility functions for logging
//

#include "tbox/Logger.h"
#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"

namespace SAMRAI {
   namespace tbox {

Logger *Logger::s_instance = NULL;

/*
 * Default Appender to print abort message and calling location to perr stream.
 */

class AbortAppender : public Logger::Appender {      

   void logMessage(const std::string &message, 
		   const std::string &filename, 
		   const int line) 
   {
      perr << "Program abort called in file ``" << filename
	   << "'' at line " << line << std::endl;
      perr << "ERROR MESSAGE: " << std::endl << message.c_str() << std::endl;
      perr << std::flush;
   }
};


/*
 * Default Appender to print a warning message and calling location to log stream.
 */
class WarningAppender : public Logger::Appender {      

   void logMessage(const std::string &message, 
		   const std::string &filename, 
		   const int line) 
   {
      plog << "Warning in file ``" << filename
	   << "'' at line " << line << std::endl;
      plog << "WARNING MESSAGE: " << std::endl << message.c_str() << std::endl;
      plog << std::flush;
   }
};


/*
 * Default Appender to print a debug message and calling location to log stream.
 */
class DebugAppender : public Logger::Appender {

   void logMessage(const std::string &message, 
		   const std::string &filename, 
		   const int line) 
   {
      plog << "Debug in file ``" << filename
	   << "'' at line " << line << std::endl;
      plog << "DEBUG MESSAGE: " << std::endl << message.c_str() << std::endl;
      plog << std::flush;
   }
};

/*
 * Default constructor for Logger singleton 
 */
Logger::Logger() : 
   d_log_warning(true),
   d_log_debug(false)
{
   /*
    * Initializers for default logging methods.
    */
   d_abort_appender = new AbortAppender();
   
   d_warning_appender = new WarningAppender();
   
   d_debug_appender = new DebugAppender();

   tbox::ShutdownRegistry::registerShutdownRoutine(free,
		   tbox::ShutdownRegistry::priorityLogger);

}

/*
 * Default destructor for Logger singleton 
 */
Logger::~Logger() 
{
   d_abort_appender.setNull();
   d_warning_appender.setNull();
   d_debug_appender.setNull();
}

void Logger::free() 
{
   if(s_instance) {
      delete s_instance;
      s_instance = static_cast<Logger *>(NULL);
   }
}


Logger *Logger::getInstance()
{
   if(s_instance == static_cast<Logger *>(NULL)) {
      s_instance = new Logger();
   }

   return s_instance;
}


/*
 * Log an abort message.
 */
void Logger::logAbort(const std::string &message, 
        			  const std::string &filename, 
        			  const int line) 
{
   d_abort_appender -> logMessage(message, filename, line);
}

/*
 * Log a warning message.
 */
void Logger::logWarning(const std::string &message, 
        			  const std::string &filename, 
        			  const int line) 
{
   if(d_log_warning) {
      d_warning_appender -> logMessage(message, filename, line);
   }
}

/*
 * Log a debug message.
 */
void Logger::logDebug(const std::string &message, 
        			  const std::string &filename, 
        			  const int line) 
{
   if(d_log_debug) {
      d_debug_appender -> logMessage(message, filename, line);
   }
}


/*
 * Set appenders.
 */
void Logger::setAbortAppender(tbox::Pointer<Appender> appender) {
   d_abort_appender = appender;
}


void Logger::setWarningAppender(tbox::Pointer<Appender> appender) {
   d_warning_appender = appender;
}


void Logger::setDebugAppender(tbox::Pointer<Appender> appender) {
   d_debug_appender = appender;
}


/*
 * Set logging state.
 */

void Logger::setWarning(bool onoff) {
   d_log_warning = onoff;
}


void Logger::setDebug(bool onoff) {
   d_log_debug = onoff;
}

}
}
