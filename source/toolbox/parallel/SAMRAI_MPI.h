//
// File:  $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/SAMRAI_MPI.h $
// Package:  SAMRAI toolbox
// Copyright:  (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:  $LastChangedRevision: 1768 $
// Modified:  $LastChangedDate: 2007-12-11 16:02:04 -0800 (Tue, 11 Dec 2007) $
// Description:  Simple utility class for interfacing with MPI
//

#ifndef included_tbox_SAMRAI_MPI
#define included_tbox_SAMRAI_MPI

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifdef HAVE_MPI
#ifndef included_mpi
#define included_mpi

/*
 * Prevent inclusion of mpi C++ bindings in mpi.h includes
 */
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif

#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#include "mpi.h"

#endif
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif


namespace SAMRAI {
   namespace tbox {


/**
 * @brief Provides C++ wrapper around MPI routines.
 *
 * Class SAMRAI_MPI groups common MPI routines into one globally-accessible
 * location.  It provides small, simple routines that are common in MPI code.
 * In some cases, the calling syntax has been simplified for convenience.
 * Moreover, there is no reason to include the preprocessor ifdef/endif
 * guards around these calls, since the MPI libraries are not called in
 * these routines if the MPI libraries are not being used (e.g., when
 * writing serial code).
 * 
 * Note that this class is a utility class to group function calls in one
 * name space (all calls are to static functions).  Thus, you should never
 * attempt to instantiate a class of type MPI; simply call the functions
 * as static functions using the MPI::function(...) syntax.
 */

struct SAMRAI_MPI
{

/**
 * MPI Types
 */
#ifdef HAVE_MPI
   typedef MPI_Comm comm;
   typedef MPI_Group group;
   typedef MPI_Request request;
   typedef MPI_Status status;
#else   
   typedef int comm;
   typedef int group;
   typedef int request;
   typedef int status;
#endif

   /**
    * MPI constants
    */
   static comm commWorld;
   static comm commNull;

   /**
    * Set boolean flag indicating whether exit or abort is called when running
    * with one processor.  Calling this function influences the behavior of
    * calls to SAMRAI_MPI::abort().  By default, the flag is true meaning that
    * abort() will be called.  Passing false means exit(-1) will be called. 
    */
   static void setCallAbortInSerialInsteadOfExit(bool flag = true); 

   /**
    * Call MPI_Abort or exit depending on whether running with one or more 
    * processes and value set by function above, if called.  The default is
    * to call exit(-1) if running with one processor and to call MPI_Abort()
    * otherwise.  This function avoids having to guard abort calls in 
    * application code.
    */
   static void abort();
   
   /**
    * Call MPI_Init.  Use of this function avoids guarding MPI init calls
    * in application code.
    */
   static void init(int* argc, char** argv[]);

   /**
    * Call MPI_Finalize.  Use of this function avoids guarding MPI finalize
    * calls in application code.
    */
   static void finalize();

   /**
    * Initialize the MPI utility class.  The MPI utility class must be
    * initialized after the call to MPI_Init or SAMRAI_MPI::init.
    */
   static void initialize();

   /**
    * Set the communicator that is used for the MPI communication routines.
    * The default communicator is MPI_COMM_WORLD.
    */
   static void setCommunicator(SAMRAI_MPI::comm communicator);

   /**
    * Get the current MPI communicator.  The default communicator is
    * MPI_COMM_WORLD.
    */
   static SAMRAI_MPI::comm getCommunicator();

   /**
    * Return the processor rank (identifier) from 0 through the number of
    * processors minus one.
    */
   static int getRank();

   /**
    * Return the number of processors (nodes).
    */
   static int getNodes();

   /**
    * Update the statistics for outgoing messages.  Statistics are
    * automatically updated for the reduction calls in MPI.
    */
   static void updateOutgoingStatistics(const int messages, const int bytes);

   /**
    * Update the statistics for incoming messages.  Statistics are
    * automatically updated for the reduction calls in MPI.
    */
   static void updateIncomingStatistics(const int messages, const int bytes);

   /**
    * Return the number of outgoing messages.
    */
   static int getOutgoingMessages();

   /**
    * Return the number of outgoing message bytes.
    */
   static int getOutgoingBytes();

   /**
    * Return the number of incoming messages.
    */
   static int getIncomingMessages();

   /**
    * Return the number of incoming message bytes.
    */
   static int getIncomingBytes();

   /**
    * Get the depth of the reduction trees given the current number
    * of MPI processors.
    */
   static int getTreeDepth();

   /**
    * Perform a global barrier across all processors.
    */
   static void barrier();

   /**
    * Perform a scalar sum reduction on a double across all nodes.  Each
    * processor contributes a value x of type double, and the sum is returned
    * from the function.
    */
   static double sumReduction(const double x);

   /**
    * Perform an array sum reduction on doubles across all nodes.  Each
    * processor contributes an array of values of type double, and the
    * element-wise sum is returned in the same array.
    */
   static void sumReduction(double *x, const int n = 1);

   /**
    * Perform a scalar sum reduction on a float across all nodes.  Each
    * processor contributes a value x of type float, and the sum is returned
    * from the function.
    */
   static float sumReduction(const float x);

   /**
    * Perform an array sum reduction on floats across all nodes.  Each
    * processor contributes an array of values of type float, and the
    * element-wise sum is returned in the same array.
    */
   static void sumReduction(float *x, const int n = 1);

   /**
    * Perform a scalar sum reduction on a dcomplex across all nodes.  Each
    * processor contributes a value x of type dcomplex, and the sum is returned
    * from the function.
    */
   static dcomplex sumReduction(const dcomplex x);

   /**
    * Perform an array sum reduction on dcomplexes across all nodes.  Each
    * processor contributes an array of values of type dcomplex, and the
    * element-wise sum is returned in the same array.
    */
   static void sumReduction(dcomplex *x, const int n = 1);

   /**
    * Perform a scalar sum reduction on an integer across all nodes.  Each
    * processor contributes a value x of type int, and the sum is returned
    * from the function.
    */
   static int sumReduction(const int x);

   /**
    * Perform an array sum reduction on integers across all nodes.
    * Each processor contributes an array of values of type int, and
    * the element-wise sum is returned in the same array.
    */
   static void sumReduction(int *x, const int n = 1);

   /**
    * Perform a scalar min reduction on a double across all nodes.  Each
    * processor contributes a value x of type double, and the minimum is
    * returned from the function.  
    *
    * If a 'rank_of_min' argument is provided, it will set it to the 
    * rank of process holding the minimum value.
    */
   static double minReduction(const double x, int *rank_of_min = NULL);

   /**
    * Perform an array min reduction on doubles across all nodes.  Each
    * processor contributes an array of values of type double, and the
    * element-wise minimum is returned in the same array.
    *
    * If a 'rank_of_min' argument is provided, it will set the array to the 
    * rank of process holding the minimum value.  Like the double argument,
    * the size of the supplied 'rank_of_min' array should be n.
    */
   static void minReduction(double *x, 
                            const int n = 1, 
                            int *rank_of_min = NULL);

   /**
    * Perform a scalar min reduction on a float across all nodes.  Each
    * processor contributes a value x of type float, and the minimum is
    * returned from the function.
    *
    * If a 'rank_of_min' argument is provided, it will set it to the 
    * rank of process holding the minimum value.
    */
   static float minReduction(const float x, int *rank_of_min = NULL);

   /**
    * Perform an array min reduction on floats across all nodes.  Each
    * processor contributes an array of values of type float, and the
    * element-wise minimum is returned in the same array.
    *
    * If a 'rank_of_min' argument is provided, it will set the array to the 
    * rank of process holding the minimum value. Like the double argument,
    * the size of the supplied 'rank_of_min' array should be n.
    */
   static void minReduction(float *x, 
                            const int n = 1, 
                            int *rank_of_min = NULL);

   /**
    * Perform a scalar min reduction on an integer across all nodes.  Each
    * processor contributes a value x of type int, and the minimum is returned
    * from the function.
    *
    * If a 'rank_of_min' argument is provided, it will set it to the 
    * rank of process holding the minimum value.
    */
   static int minReduction(const int x, int *rank_of_min = NULL);

   /**
    * Perform an array min reduction on integers across all nodes.
    * Each processor contributes an array of values of type int, and
    * the element-wise minimum is returned in the same array.
    *
    * If a 'rank_of_min' argument is provided, it will set the array to the 
    * rank of process holding the minimum value. Like the double argument,
    * the size of the supplied 'rank_of_min' array should be n.
    */
   static void minReduction(int *x, 
                            const int n = 1, 
                            int *rank_of_min = NULL);

   /**
    * Perform a scalar max reduction on a double across all nodes.  Each
    * processor contributes a value x of type double, and the maximum is
    * returned from the function.
    *
    * If a 'rank_of_max' argument is provided, it will set it to the 
    * rank of process holding the maximum value.
    */
   static double maxReduction(const double x, int *rank_of_max = NULL);

   /**
    * Perform an array max reduction on doubles across all nodes.  Each
    * processor contributes an array of values of type double, and the
    * element-wise maximum is returned in the same array.
    *
    * If a 'rank_of_max' argument is provided, it will set the array to the 
    * rank of process holding the maximum value. Like the double argument,
    * the size of the supplied 'rank_of_max' array should be n.
    */
   static void maxReduction(double *x, 
                            const int n = 1, 
                            int *rank_of_max = NULL);

   /**
    * Perform a scalar max reduction on a float across all nodes.  Each
    * processor contributes a value x of type float, and the maximum is
    * returned from the function.
    *
    * If a 'rank_of_max' argument is provided, it will set it to the 
    * rank of process holding the maximum value.
    */
   static float maxReduction(const float x, int *rank_of_max = NULL);

   /**
    * Perform an array max reduction on floats across all nodes.  Each
    * processor contributes an array of values of type float, and the
    * element-wise maximum is returned in the same array.
    *
    * If a 'rank_of_max' argument is provided, it will set the array to the 
    * rank of process holding the maximum value. Like the double argument,
    * the size of the supplied 'rank_of_max' array should be n.
    */
   static void maxReduction(float *x, 
                            const int n = 1, 
                            int *rank_of_max = NULL);

   /**
    * Perform a scalar max reduction on an integer across all nodes.  Each
    * processor contributes a value x of type int, and the maximum is returned
    * from the function.
    *
    * If a 'rank_of_max' argument is provided, it will set it to the 
    * rank of process holding the maximum value.
    */
   static int maxReduction(const int x, int *rank_of_max = NULL);

   /**
    * Perform an array max reduction on integers across all nodes.
    * Each processor contributes an array of values of type int, and
    * the element-wise maximum is returned in the same array.
    *
    * If a 'rank_of_max' argument is provided, it will set the array to the 
    * rank of process holding the maximum value. Like the double argument,
    * the size of the supplied 'rank_of_max' array should be n.
    */
   static void maxReduction(int *x, 
                            const int n = 1, 
                            int *rank_of_max = NULL);

   /**
    * Perform an all-to-one sum reduction on an integer array.
    * The final result is only available on the root processor.
    */
   static void allToOneSumReduction(int *x, const int n, const int root = 0);

   /**
    * Broadcast integer from specified root process to all other processes.
    * All processes other than root, receive a copy of the integer value.
    */
   static int bcast(const int x, const int root);

   /**
    * Broadcast integer array from specified root processor to all other
    * processors.  For the root processor, "array" and "length"
    * are treated as const.
    */
   static void bcast(int *x, int &length, const int root);

   /**
    * Broadcast char array from specified root processor to all other
    * processors.  For the root processor, "array" and "length"
    * are treated as const.
    */
   static void bcast(char *x, int &length, const int root);

   /*!
    * @brief This function sends an MPI message with an integer 
    * array to another processer.
    *
    * If the receiving processor knows in advance the length 
    * of the array, use "send_length = false;"  otherwise, 
    * this processor will first send the length of the array, 
    * then send the data.  This call must be paired  with a 
    * matching call to SAMRAI_MPI::recv.
    *
    * @param buf Pointer to integer array buffer with length integers.
    * @param length Number of integers in buf that we want to send.
    * @param receiving_proc_number Receiving processor number.
    * @param send_length Optional boolean argument specifiying if 
    * we first need to send a message with the array size.
    * Default value is true.
    * @param tag Optional integer argument specifying an integer tag
    * to be sent with this message.  Default tag is 0.
    */

   static void send(const int *buf, 
                    const int length, 
                    const int receiving_proc_number, 
                    const bool send_length = true, 
                    int tag = -1);

   /*!
    * @brief This function sends an MPI message with an array of bytes
    * (MPI_BYTES) to receiving_proc_number.
    *
    * This call must be paired with a matching call to SAMRAI_MPI::recvBytes.
    *
    * @param buf Void pointer to an array of number_bytes bytes to send.
    * @param number_bytes Integer number of bytes to send.
    * @param receiving_proc_number Receiving processor number.
    */
  static void sendBytes(const void *buf,
                        const int number_bytes, 
                        const int receiving_proc_number);

   /*!
    * @brief This function receives an MPI message with an array of
    * max size number_bytes (MPI_BYTES) from any processer.
    *
    * This call must be paired with a matching call to SAMRAI_MPI::sendBytes.
    *               
    * This function returns the processor number of the sender.
    *
    * @param buf Void pointer to a buffer of size number_bytes bytes.
    * @param number_bytes Integer number specifing size of buf in bytes.
    */
   static int recvBytes(void *buf, int number_bytes);

   /*!
    * @brief This function receives an MPI message with an integer 
    * array from another processer.
    *
    * If this processor knows in advance the length of the array,
    * use "get_length = false;" otherwise, the sending processor 
    * will first send the length of the array, then send the data.
    * This call must be paired with a matching call to SAMRAI_MPI::send.
    *
    * @param buf Pointer to integer array buffer with capacity of
    * length integers.
    * @param length Maximum number of integers that can be stored in
    * buf.
    * @param sending_proc_number Processor number of sender.
    * @param get_length Optional boolean argument specifiying if 
    * we first need to send a message to determine the array size.
    * Default value is true.
    * @param tag Optional integer argument specifying a tag which
    * must be matched by the tag of the incoming message. Default
    * tag is 0.
    */
   static void recv(int *buf, 
                    int &length, 
                    const int sending_proc_number,
                    const bool get_length = true, 
                    int tag = -1);
//@}
   /**
    * Each processor sends an array of integers or doubles to all other
    * processors; each processor's array may differ in length.
    * The x_out array must be pre-allocated to the correct length
    * (this is a bit cumbersome, but is necessary to avoid th allGather
    * function from allocating memory that is freed elsewhere).
    * To properly preallocate memory, before calling this method, call
    *
    *   size_out = SAMRAI_MPI::sumReduction(size_in)
    *
    * then allocate the x_out array.
    */
   static void allGather(const int *x_in, int size_in, 
                         int *x_out, int size_out);
   static void allGather(const double *x_in, int size_in, 
                         double *x_out, int size_out);
//@}



//@{
   /**
    * Each processor sends every other processor an integer or double.
    * The x_out array should be preallocated to a length equal
    * to the number of processors.
    */
   static void allGather(int x_in, int *x_out);
   static void allGather(double x_in, double *x_out);

//@}


private:
   static SAMRAI_MPI::comm     s_communicator;
   static int      s_outgoing_messages;
   static int      s_outgoing_bytes;
   static int      s_incoming_messages;
   static int      s_incoming_bytes;
   static int      s_initialized;

   static bool     s_call_abort_in_serial_instead_of_exit;

   /**
    * Performs common functions needed by some of the allToAll methods.
    */
   static void allGatherSetup(int size_in, int size_out, 
                              int *&rcounts, int *&disps);

   //@{
   //@name Structs for passing arguments to MPI
   struct DoubleIntStruct { double d; int i; };
   struct FloatIntStruct { float f; int i; };
   struct IntIntStruct { int j; int i; };
   //@}

};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/SAMRAI_MPI.I"
#endif
#endif
