/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for logging
 *
 ************************************************************************/

#include "SAMRAI/tbox/Logger.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace tbox {

Logger * Logger::s_instance = NULL;

StartupShutdownManager::Handler
Logger::s_finalize_handler(
   0,
   0,
   0,
   Logger::finalizeCallback,
   tbox::StartupShutdownManager::priorityLogger);

/*
 * Default Appender to print abort message and calling location to perr stream.
 */

class AbortAppender:public Logger::Appender
{

   void logMessage(
      const std::string& message,
      const std::string& filename,
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
class WarningAppender:public Logger::Appender
{

   void logMessage(
      const std::string& message,
      const std::string& filename,
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
class DebugAppender:public Logger::Appender
{

   void logMessage(
      const std::string& message,
      const std::string& filename,
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
Logger::Logger():
   d_log_warning(true),
   d_log_debug(false)
{
   /*
    * Initializers for default logging methods.
    */
   d_abort_appender = new AbortAppender();

   d_warning_appender = new WarningAppender();

   d_debug_appender = new DebugAppender();

}

/*
 * Default destructor for Logger singleton
 */
Logger::~Logger()
{
   d_abort_appender.reset();
   d_warning_appender.reset();
   d_debug_appender.reset();
}

void Logger::finalizeCallback()
{
   if (s_instance) {
      delete s_instance;
      s_instance = static_cast<Logger *>(NULL);
   }
}

Logger *Logger::getInstance()
{
   if (s_instance == static_cast<Logger *>(NULL)) {
      s_instance = new Logger();
   }

   return s_instance;
}

/*
 * Log an abort message.
 */
void Logger::logAbort(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   d_abort_appender->logMessage(message, filename, line);
}

/*
 * Log a warning message.
 */
void Logger::logWarning(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   if (d_log_warning) {
      d_warning_appender->logMessage(message, filename, line);
   }
}

/*
 * Log a debug message.
 */
void Logger::logDebug(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   if (d_log_debug) {
      d_debug_appender->logMessage(message, filename, line);
   }
}

/*
 * Set appenders.
 */
void Logger::setAbortAppender(
   tbox::Pointer<Appender> appender) {
   d_abort_appender = appender;
}

void Logger::setWarningAppender(
   tbox::Pointer<Appender> appender) {
   d_warning_appender = appender;
}

void Logger::setDebugAppender(
   tbox::Pointer<Appender> appender) {
   d_debug_appender = appender;
}

/*
 * Set logging state.
 */

void Logger::setWarning(
   bool onoff) {
   d_log_warning = onoff;
}

void Logger::setDebug(
   bool onoff) {
   d_log_debug = onoff;
}

}
}
