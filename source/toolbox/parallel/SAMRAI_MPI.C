//
// File:  $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/SAMRAI_MPI.C $
// Package:  SAMRAI toolbox
// Copyright:  (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:  $LastChangedRevision: 1704 $
// Modified:  $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:  Simple utility class for interfacing with MPI
//


#include <stdlib.h>

#include <string>
using namespace std;

#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/SAMRAI_MPI.I"
#endif

#ifdef __INSURE__
/*
 * These are defined in mpich mpi.h and break the insure compile.
 * This may impact Globus in some way, at least from the comments
 * in the mpi.h header file.  Why mpich externs something that is 
 * not defined in the mpich is confusing and probably just broken.
 */
int MPICHX_TOPOLOGY_DEPTHS;
int MPICHX_TOPOLOGY_COLORS;
int MPICHX_PARALLELSOCKETS_PARAMETERS;
#endif



namespace SAMRAI {
   namespace tbox {

SAMRAI_MPI::comm SAMRAI_MPI::s_communicator      = (SAMRAI_MPI::comm) 0;
int      SAMRAI_MPI::s_outgoing_messages = 0;
int      SAMRAI_MPI::s_outgoing_bytes    = 0;
int      SAMRAI_MPI::s_incoming_messages = 0;
int      SAMRAI_MPI::s_incoming_bytes    = 0;
int      SAMRAI_MPI::s_initialized       = 0; 

#if HAVE_MPI
SAMRAI_MPI::comm SAMRAI_MPI::commWorld = MPI_COMM_WORLD;
SAMRAI_MPI::comm SAMRAI_MPI::commNull = MPI_COMM_NULL;
#else
SAMRAI_MPI::comm SAMRAI_MPI::commWorld = 0;
SAMRAI_MPI::comm SAMRAI_MPI::commNull = -1;
#endif

bool SAMRAI_MPI::s_call_abort_in_serial_instead_of_exit = true;

/*
**************************************************************************
*                                                                        *
* Abort the program.                                                     *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(bool flag)
{
   s_call_abort_in_serial_instead_of_exit = flag;
}

void SAMRAI_MPI::abort()
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      MPI_Abort(s_communicator, -1);
   } else {
      if (s_call_abort_in_serial_instead_of_exit) {
         ::abort();
      } else {
         exit(-1);
      }
   }
#else
   if (s_call_abort_in_serial_instead_of_exit) {
      ::abort();
   } else {
      exit(-1);
   }
#endif

}

/*
**************************************************************************
*                                                                        *
* Initialize the static data in the MPI utility class.  This must be     *
* called after MPI_Init to ensure that the MPI_COMM_WORLD structure      *
* has been initialized.                                                  *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::initialize()
{
   if (!s_initialized) {
      s_communicator      = SAMRAI_MPI::commWorld;
      s_outgoing_messages = 0;
      s_outgoing_bytes    = 0;
      s_incoming_messages = 0;
      s_incoming_bytes    = 0;
      s_initialized = 1;
   }
}

/*
**************************************************************************
*                                                                        *
* Tree depth calculation for tracking the * number of message sends      *
* and receives and the number of bytes.                                  *
*                                                                        *
**************************************************************************
*/

int SAMRAI_MPI::getTreeDepth()
{
   int depth = 0;
   const int nnodes = getNodes();
   while ((1 << depth) < nnodes) {
      depth++;
   }
   return(depth);
}

/*
**************************************************************************
*                                                                        *
* Perform a global barrier across all processors.                        *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::barrier()
{
#ifdef HAVE_MPI
   (void) MPI_Barrier(s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, 0);
   updateIncomingStatistics(tree, 0);
#endif
}
 
/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double SAMRAI_MPI::sumReduction(const double x)
{
   double recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::sumReduction(double *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double *send = new double[n];
      memcpy(send, x, n*sizeof(double));
      MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float SAMRAI_MPI::sumReduction(const float x)
{
   float recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      float send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_FLOAT, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::sumReduction(float *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      float *send = new float[n];
      memcpy(send, x, n*sizeof(float));
      MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar dcomplex.                                   *
*                                                                        *
**************************************************************************
*/

dcomplex SAMRAI_MPI::sumReduction(const dcomplex x)
{
   dcomplex recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double xreal[2];
      double send[2];
      send[0] = xreal[0] = real(x); 
      send[1] = xreal[1] = imag(x);
      MPI_Allreduce(send, xreal, 2, MPI_DOUBLE, MPI_SUM, s_communicator);
      recv = dcomplex(xreal[0], xreal[1]);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, 2*sizeof(double));
      updateIncomingStatistics(tree, 2*sizeof(double));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of dcomplex.                                *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::sumReduction(dcomplex *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int nrvals = 2*n;
      double *xreal = new double[nrvals];
      double *send = new double[nrvals];
      for (int i = 0; i < n; i++) {
         xreal[2*i]   = real(x[i]); 
         xreal[2*i+1] = imag(x[i]); 
      }
      memcpy(send, xreal, nrvals*sizeof(double));
      MPI_Allreduce(send, xreal, nrvals, MPI_DOUBLE, MPI_SUM, s_communicator);
      for (int j = 0; j < n; j++) {
         x[j] = dcomplex(xreal[2*j], xreal[2*j+1]);
      }
      delete [] send;
      delete [] xreal;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int SAMRAI_MPI::sumReduction(const int x)
{
   int recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of integers.                                z
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::sumReduction(int *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int *send = new int[n];
      memcpy(send, x, n*sizeof(int));
      MPI_Allreduce(send, x, n, MPI_INT, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double SAMRAI_MPI::minReduction(const double x, int *rank_of_min)
{
   double rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         double send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MIN, s_communicator);
      } else {
         DoubleIntStruct recv;
         DoubleIntStruct send;
         send.d = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_DOUBLE_INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.d;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::minReduction(double *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         double *send = new double[n];
         memcpy(send, x, n*sizeof(double));
         MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_DOUBLE_INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            rank_of_min[i] = send[i].i;
         }
         delete [] recv;
         delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float SAMRAI_MPI::minReduction(const float x, int *rank_of_min)
{
   float rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         float send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MIN, s_communicator);
      }
      else {
         FloatIntStruct recv;
         FloatIntStruct send;
         send.f = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_FLOAT_INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.f;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::minReduction(float *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         float *send = new float[n];
         memcpy(send, x, n*sizeof(float));
         MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_FLOAT_INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            rank_of_min[i] = send[i].i;
         }
         delete [] recv;
         delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int SAMRAI_MPI::minReduction(const int x, int *rank_of_min)
{
   int rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         int send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MIN, s_communicator);
      }
      else {
         IntIntStruct recv;
         IntIntStruct send;
         send.j = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_2INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.j;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of integers.                                *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::minReduction(int *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         int *send = new int[n];
         memcpy(send, x, n*sizeof(int));
         MPI_Allreduce(send, x, n, MPI_INT, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_2INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            rank_of_min[i] = send[i].i;
         }
         delete [] recv;
         delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);

#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double SAMRAI_MPI::maxReduction(const double x, int *rank_of_max)
{
   double rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         double send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MAX, s_communicator);
      } else {
         DoubleIntStruct recv;
         DoubleIntStruct send;
         send.d = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_DOUBLE_INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.d;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::maxReduction(double *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         double *send = new double[n];
         memcpy(send, x, n*sizeof(double));
         MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_DOUBLE_INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            rank_of_max[i] = send[i].i;
         }
         delete [] recv;
         delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float SAMRAI_MPI::maxReduction(const float x, int *rank_of_max)
{
   float rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         float send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MAX, s_communicator);
      }
      else {
         FloatIntStruct recv;
         FloatIntStruct send;
         send.f = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_FLOAT_INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.f;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::maxReduction(float *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         float *send = new float[n];
         memcpy(send, x, n*sizeof(float));
         MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_FLOAT_INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            rank_of_max[i] = send[i].i;
         }
         delete [] recv;
         delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int SAMRAI_MPI::maxReduction(const int x, int *rank_of_max)
{
   int rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         int send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MAX, s_communicator);
      }
      else {
         IntIntStruct recv;
         IntIntStruct send;
         send.j = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_2INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.j;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of integers.                                *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::maxReduction(int *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         int *send = new int[n];
         memcpy(send, x, n*sizeof(int));
         MPI_Allreduce(send, x, n, MPI_INT, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_2INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            rank_of_max[i] = send[i].i;
         }
	 delete [] recv;
	 delete [] send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-one sum reduction on integer array.                             *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::allToOneSumReduction(int *x, const int n, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int *send = new int[n];
      memcpy(send, x, n*sizeof(int));
      MPI_Reduce(send, x, n, MPI_INT, MPI_SUM, root, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();

      /*
       * note: following probably isn't correct, since some processors
       * both send and receive, and some only send; this is of course
       * dependent on the MPI implementation and/or the hardware.
       */
      if (getRank() == root) {
         updateOutgoingStatistics(tree, n*sizeof(int));
      } else {
         updateIncomingStatistics(tree, n*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}


/*
**************************************************************************
*                                                                        *
* Broadcast for scalar integer.                                          *
*                                                                        *
**************************************************************************
*/

int SAMRAI_MPI::bcast(const int x, const int root)
{
   int recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      (void) MPI_Bcast(&recv, 1, MPI_INT, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, sizeof(int));
      } else {
         updateIncomingStatistics(tree, sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Broadcast for integer array from root processor to all other procs.    *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::bcast(int *x, int &length, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {

      (void) MPI_Bcast((void*)x, length, MPI_INT, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, length*sizeof(int));
      } else {
         updateIncomingStatistics(tree, length*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
}

/*
**************************************************************************
*                                                                        *
* Broadcast for char array from root processor to all other processors.  *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::bcast(char *x, int &length, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {

      (void) MPI_Bcast((void*)x, length, MPI_BYTE, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, length*sizeof(int));
      } else {
         updateIncomingStatistics(tree, length*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
}

/*
**************************************************************************
*                                                                        *
* Send integer array to another processor.                               *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::send(const int *buf, 
                    const int length, 
                    const int receiving_proc_number, 
                    const bool send_length, 
                    int tag)
{
#ifdef HAVE_MPI
   tag = (tag >= 0) ? tag : 0;
   int size = length;
   if (send_length) {
      MPI_Send(&size, 1, MPI_INT, receiving_proc_number, tag, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
   }
   MPI_Send((void*)buf, 
            length, 
            MPI_INT, 
            receiving_proc_number, 
            tag, 
            s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, length * sizeof(int));
#endif
}

/*
**************************************************************************
*                                                                        *
* Receive integer array from another processor.                          *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::recv(int *buf, 
                    int &length, 
                    const int sending_proc_number,
		    const bool get_length, 
                    int tag)
{
#ifdef HAVE_MPI
   MPI_Status status;
   tag = (tag >= 0) ? tag : 0;
   if (get_length) {
      MPI_Recv(&length, 
               1, 
               MPI_INT, 
               sending_proc_number, 
               tag, 
               s_communicator,
               &status);
      const int tree = getTreeDepth();
      updateIncomingStatistics(tree, sizeof(int));
   }
   MPI_Recv((void*)buf, 
            length, 
            MPI_INT, 
            sending_proc_number, 
            tag, 
            s_communicator, 
            &status);
   const int tree = getTreeDepth();
   updateIncomingStatistics(tree, length*sizeof(int));
#endif
}


/*
*************************************************************************
*									*
* Send an array of number_bytes bytes from this processer to            *
* receiving_proc.                                                       *
* This call must be paired with a matching call to SAMRAI_MPI::recvBytes. *
*									*
*************************************************************************
*/

void SAMRAI_MPI::sendBytes(const void *buf, 
                         const int number_bytes, 
                         const int receiving_proc_number)
{
#ifdef HAVE_MPI

   MPI_Send((void*)buf, 
            number_bytes, 
            MPI_BYTE, 
            receiving_proc_number, 
            0, 
            s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, number_bytes * sizeof(char));
#endif
}



/*
*************************************************************************
*									*
* Receive an array of bytes of max size number_bytes bytes from any     *
* processer.                                                            *
* This call must be paired with a matching call to SAMRAI_MPI::sendBytes. *
*                                                                       *
* Returns the processor number of the sender.                           *
*									*
*************************************************************************
*/

int SAMRAI_MPI::recvBytes(void *buf, 
                        int number_bytes) 
{
   int rval = 0;
#ifdef HAVE_MPI
   MPI_Status status;
   MPI_Recv(buf, 
            number_bytes, 
            MPI_BYTE, 
            MPI_ANY_SOURCE, 
            MPI_ANY_TAG, 
            s_communicator, 
            &status);

   const int tree = getTreeDepth();
   updateIncomingStatistics(tree, number_bytes * sizeof(char));
   rval = status.MPI_SOURCE;
#endif

   return rval;
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of arrays of integers; each processor's            *
* array can be of a different length                                     *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::allGather(
   const int *x_in, int size_in, int *x_out, int size_out)
{
#ifdef HAVE_MPI
   int* rcounts = (int*)NULL;
   int* disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   MPI_Allgatherv((void*)x_in, size_in, MPI_INT,
                  x_out, rcounts, disps, MPI_INT, s_communicator);

   if (rcounts) {
      delete [] rcounts;
   }
   if (disps) {
      delete [] disps;
   }
#else
   NULL_USE(x_in);
   NULL_USE(size_in);
   NULL_USE(x_out);
   NULL_USE(size_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* all-to-all exchange of arrays of doubles; each processor's             *
* array can be of a different length                                     *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::allGather(
   const double *x_in, int size_in, double *x_out, int size_out)
{
#ifdef HAVE_MPI
   int *rcounts = (int*)NULL;
   int *disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   MPI_Allgatherv((void*)x_in, size_in, MPI_DOUBLE,
                  x_out, rcounts, disps, MPI_DOUBLE, s_communicator);

   if (rcounts) {
      delete [] rcounts;
   }
   if (disps) {
      delete [] disps;
   }
#else
   NULL_USE(x_in);
   NULL_USE(size_in);
   NULL_USE(x_out);
   NULL_USE(size_out);
#endif
}

/*
*************************************************************************
*                                                                       *
* common setup funtion for all-to-all functions                         *
*                                                                       *
*************************************************************************
*/

void SAMRAI_MPI::allGatherSetup(
   int size_in, int size_out, int *&rcounts, int *&disps)
{
#ifdef HAVE_MPI
   int np = getNodes();
   rcounts = new int[np];
   disps = new int[np];

   /* figure out where where each processor's input will be placed */
   allGather(size_in, rcounts);

   disps[0] = 0;
   for (int p = 1; p < np; ++p) {
      disps[p] = disps[p-1] + rcounts[p-1];
   }

   /* verify that the x_out array is the appropriate size! */
   int c = 0;
   for (int x = 0; x < np; ++x) {
      c += rcounts[x];
   }
   if (c != size_out) {
      TBOX_ERROR("SAMRAI_MPI::allGatherSetup error..." 
                 << "\n   size_out =" << size_out << "appears to be incorrect; "
                 << "should be: " << c << endl);
   }
#else
   NULL_USE(size_in);
   NULL_USE(size_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of a double.                                       *
*                                                                        *
**************************************************************************
*/

void SAMRAI_MPI::allGather(double x_in, double *x_out)
{
#ifdef HAVE_MPI
   MPI_Allgather(&x_in, 1, MPI_DOUBLE, x_out, 1, MPI_DOUBLE, s_communicator);
#else
   NULL_USE(x_in);
   NULL_USE(x_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of an integer.                                     *
*                                                                        *
* ************************************************************************
*/

void SAMRAI_MPI::allGather(int x_in, int *x_out)
{
#ifdef HAVE_MPI
   MPI_Allgather(&x_in, 1, MPI_INT, x_out, 1, MPI_INT, s_communicator);
#else
   NULL_USE(x_in);
   NULL_USE(x_out);
#endif
}

}
}
