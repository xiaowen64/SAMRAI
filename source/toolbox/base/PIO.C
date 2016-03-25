//
// File:	PIO.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Parallel I/O classes pout, perr, and plog and control class
//

#include "tbox/PIO.h"
#include <stdio.h>
#include "tbox/MPI.h"
#include "tbox/ParallelBuffer.h"

#ifndef NULL
#define NULL 0
#endif

namespace SAMRAI {
   namespace tbox {

int       PIO::s_rank       = -1;
ofstream* PIO::s_filestream = NULL;

/*
*************************************************************************
*									*
* Define the parallel buffers and the associated ostream objects.	*
*									*
*************************************************************************
*/

static ParallelBuffer pout_buffer;
static ParallelBuffer perr_buffer;
static ParallelBuffer plog_buffer;

ostream pout(&pout_buffer);
ostream perr(&perr_buffer);
ostream plog(&plog_buffer);

/*
*************************************************************************
*									*
* Initialie the parallel I/O streams.  This routine must be called	*
* before pout, perr, and plog are used for output but after MPI has	*
* been initialized.  By default, logging is disabled.			*
*									*
*************************************************************************
*/

void PIO::initialize()
{
   s_rank       = MPI::getRank();
   s_filestream = NULL;
   
   /*
    * Initialize the standard parallel output stream
    */

   pout_buffer.setActive(s_rank == 0);
   pout_buffer.setPrefixString(string());
   pout_buffer.setOutputStream1(&cout);
   pout_buffer.setOutputStream2(NULL);

   /*
    * Initialize the error parallel output stream
    */

   char buffer[16];
   sprintf(buffer, "P=%05d:", s_rank);

   perr_buffer.setActive(true);
   perr_buffer.setPrefixString(buffer);
   perr_buffer.setOutputStream1(&cerr);
   perr_buffer.setOutputStream2(NULL);

   /*
    * Initialize the parallel log file (disabled by default)
    */

   plog_buffer.setActive(false);
   plog_buffer.setPrefixString(string());
   plog_buffer.setOutputStream1(NULL);
   plog_buffer.setOutputStream2(NULL);
}

/*
*************************************************************************
*									*
* Close the output streams.  Flush both cout and cerr.  If logging,	*
* then flush and close the log stream.					*
*									*
*************************************************************************
*/

void PIO::finalize()
{
   cout.flush();
   cerr.flush();
   shutdownFilestream();
}

/*
*************************************************************************
*									*
* If the log file stream is open, then shut down the filestream.  Close	*
* and flush the channel and disconnect the output stream buffers.	*
*									*
*************************************************************************
*/

void PIO::shutdownFilestream()
{
   if (s_filestream) {
      s_filestream->flush();
      s_filestream->close();

      delete s_filestream;
      s_filestream = NULL;

      pout_buffer.setOutputStream2(NULL);
      perr_buffer.setOutputStream2(NULL);
      plog_buffer.setOutputStream1(NULL);
      plog_buffer.setActive(false);
   }
}

/*
*************************************************************************
*									*
* Log messages for node zero only.  If a log stream was open, close	*
* it.  If this is node zero, then open a new log stream and set the	*
* appropriate buffer streams to point to the log file.			*
*									*
*************************************************************************
*/

void PIO::logOnlyNodeZero(const string &filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * If this is node zero, then open the log stream and redirect output
    */

   if (s_rank == 0) {
      s_filestream = new ofstream(filename.c_str());
      if (!(*s_filestream)) {
         delete s_filestream;
         s_filestream = NULL;
         perr << "PIO: Could not open log file ``" << filename.c_str() << "''\n";
      } else {
         pout_buffer.setOutputStream2(s_filestream);
         perr_buffer.setOutputStream2(s_filestream);
         plog_buffer.setOutputStream1(s_filestream);
         plog_buffer.setActive(true);
      }
   }
}

/*
*************************************************************************
*									*
* Log messages for all nodes.  If a log stream was open, the close it.	*
* Open a log stream on every processor.  The filename for the log file	*
* will be appended with the processor number.				*
*									*
*************************************************************************
*/

void PIO::logAllNodes(const string &filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * Open the log stream and redirect output
    */

   char *buffer = new char[filename.length() + 16];
   sprintf(buffer, "%s.%05d", filename.c_str(), s_rank);
   s_filestream = new ofstream(buffer);

   if (!(*s_filestream)) {
      delete s_filestream;
      s_filestream = NULL;
      perr << "PIO: Could not open log file ``" << buffer << "''\n";
   } else {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }

   delete [] buffer;
}

/*
*************************************************************************
*									*
* Suspend logging of data to the file stream.  This does not close the	*
* filestream (assuming it is open) but just disables logging.		*
*									*
*************************************************************************
*/

void PIO::suspendLogging()
{
   pout_buffer.setOutputStream2(NULL);
   perr_buffer.setOutputStream2(NULL);
   plog_buffer.setOutputStream1(NULL);
   plog_buffer.setActive(false);
}

/*
*************************************************************************
*									*
* Resume logging of the file stream (assuming it was open).  If the	*
* file stream is NULL, then do nothing.					*
*									*
*************************************************************************
*/

void PIO::resumeLogging()
{
   if (s_filestream) {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}


}
}
