/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility class for logging.
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 * Log an abort message.
 */
SAMRAI_INLINE_KEYWORD
void
Logger::logAbort(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   d_abort_appender->logMessage(message, filename, line);
}

/*
 * Log a warning message.
 */
SAMRAI_INLINE_KEYWORD
void
Logger::logWarning(
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
SAMRAI_INLINE_KEYWORD
void
Logger::logDebug(
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
SAMRAI_INLINE_KEYWORD
void
Logger::setAbortAppender(
   const boost::shared_ptr<Appender>& appender) {
   d_abort_appender = appender;
}

SAMRAI_INLINE_KEYWORD
void
Logger::setWarningAppender(
   const boost::shared_ptr<Appender>& appender) {
   d_warning_appender = appender;
}

SAMRAI_INLINE_KEYWORD
void
Logger::setDebugAppender(
   const boost::shared_ptr<Appender>& appender) {
   d_debug_appender = appender;
}

/*
 * Set logging state.
 */

SAMRAI_INLINE_KEYWORD
void
Logger::setWarning(
   bool onoff) {
   d_log_warning = onoff;
}

SAMRAI_INLINE_KEYWORD
void
Logger::setDebug(
   bool onoff) {
   d_log_debug = onoff;
}

}
}
