//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/parallel/SAMRAI_MPI.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Simple utility class for interfacing with MPI
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {



inline
void SAMRAI_MPI::init(int* argc, char** argv[]) 
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
void SAMRAI_MPI::finalize()
{
#ifdef HAVE_MPI
   MPI_Finalize();
#endif
}

inline
void SAMRAI_MPI::setCommunicator(SAMRAI_MPI::comm communicator)
{
   s_communicator = communicator;
}

inline
SAMRAI_MPI::comm SAMRAI_MPI::getCommunicator()
{
   return(s_communicator);
}

inline
int SAMRAI_MPI::getRank()
{
   int myid = 0;
#ifdef HAVE_MPI
   MPI_Comm_rank(s_communicator, &myid);
#endif
   return(myid);
}

inline
int SAMRAI_MPI::getNodes()
{
   int nodes = 1;
#ifdef HAVE_MPI
   MPI_Comm_size(s_communicator, &nodes);
#endif
   return(nodes);
}

inline
void SAMRAI_MPI::updateOutgoingStatistics(const int messages, const int bytes)
{
   s_outgoing_messages += messages;
   s_outgoing_bytes    += bytes;
}

inline
void SAMRAI_MPI::updateIncomingStatistics(const int messages, const int bytes)
{
   s_incoming_messages += messages;
   s_incoming_bytes    += bytes;
}

inline
int SAMRAI_MPI::getOutgoingMessages()
{
   return(s_outgoing_messages);
}

inline
int SAMRAI_MPI::getOutgoingBytes()
{
   return(s_outgoing_bytes);
}

inline
int SAMRAI_MPI::getIncomingMessages()
{
   return(s_incoming_messages);
}

inline
int SAMRAI_MPI::getIncomingBytes()
{
   return(s_incoming_bytes);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
