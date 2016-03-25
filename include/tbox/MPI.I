//
// File:	MPI.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Simple utility class for interfacing with MPI
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {



inline
void MPI::init(int* argc, char** argv[]) 
{
#ifdef HAVE_MPI
   int test;
   // Determine if MPI has been initialized and init if needed
   MPI_Initialized(&test);
   if(!test)
      MPI_Init(argc,argv);
#else
   (void) argc;
   (void) argv;
#endif
}

inline
void MPI::finalize()
{
#ifdef HAVE_MPI
   MPI_Finalize();
#endif
}

inline
void MPI::setCommunicator(MPI::comm communicator)
{
   s_communicator = communicator;
}

inline
MPI::comm MPI::getCommunicator()
{
   return(s_communicator);
}

inline
int MPI::getRank()
{
   int myid = 0;
#ifdef HAVE_MPI
   MPI_Comm_rank(s_communicator, &myid);
#endif
   return(myid);
}

inline
int MPI::getNodes()
{
   int nodes = 1;
#ifdef HAVE_MPI
   MPI_Comm_size(s_communicator, &nodes);
#endif
   return(nodes);
}

inline
void MPI::updateOutgoingStatistics(const int messages, const int bytes)
{
   s_outgoing_messages += messages;
   s_outgoing_bytes    += bytes;
}

inline
void MPI::updateIncomingStatistics(const int messages, const int bytes)
{
   s_incoming_messages += messages;
   s_incoming_bytes    += bytes;
}

inline
int MPI::getOutgoingMessages()
{
   return(s_outgoing_messages);
}

inline
int MPI::getOutgoingBytes()
{
   return(s_outgoing_bytes);
}

inline
int MPI::getIncomingMessages()
{
   return(s_incoming_messages);
}

inline
int MPI::getIncomingBytes()
{
   return(s_incoming_bytes);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
