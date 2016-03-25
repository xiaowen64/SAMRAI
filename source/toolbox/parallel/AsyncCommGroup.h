/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/AsyncCommGroup.h $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: All-to-one and one-to-all communication using a tree.
 */

#ifndef included_tbox_AsyncCommGroup
#define included_tbox_AsyncCommGroup

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif

#ifndef included_tbox_RelaunchableJob
#include "tbox/RelaunchableJob.h"
#endif

#ifndef included_tbox_SAMRAI_MPI
#include "tbox/SAMRAI_MPI.h"
#endif

#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif

namespace SAMRAI {
   namespace tbox {

/*!
 * @brief Supports many-to-one and one-to-many asynchronous
 * communication operations within a given group of processes
 * by sending messages along the branches of a conceptual tree.
 *
 * This class was created to perform certain group communications
 * without using MPI global communications, which require
 * creating new MPI communicators (can be expensive) and
 * does not support asynchronous operations.
 *
 * The supported communications are asynchronous in that
 * you can start one and wait for it or check back on it
 * occassionally until it completes.
 * Asynchronous operations in conjunction with other groups
 * can be done by using a AsyncCommStage to allocate
 * the groups and to check for completed communications.
 *
 * Supported operations are currently broadcast, gather and
 * sum reduce.  Only integer data is supported.
 *
 * A tree is an acyclic graph in which a node at position pos
 * has nchild children, and the following positions
 * for its
 * - parent: (pos-1)/nchild
 * - first (oldest) child: pos*nchild+1
 * - last (youngest) child: (pos+1)*nchild
 *
 * For example, nchild=2 corresponds to a binary tree.
 *
 * Communication is done by sending messages toward the
 * root (for all-to-one operations) or leaves (for
 * one-to-all operations).  For the former, we receive
 * data from the children and send to the parent.  For
 * the latter, we receive from the parent and send to
 * the children.  Thus every communication involves a
 * receive and a send (except at the root and leaf
 * nodes of the tree).
 *
 * Using a tree generally gives
 * better performance than having all processes in the
 * the tree communicate directly with the root process.
 * Using MPI communicators corresponding to the groups
 * may faster than using this class, but the cost of
 * creating MPI communicators MAY be expensive.
 *
 * This class supports communication and uses MPI for
 * message passing.  If MPI is disabled, the job of this
 * class disappears and the class is effectively empty.  
 * The public interfaces still remain so the class
 * can compile, but the implementations are trivial.
 */
class AsyncCommGroup
{

public:


   /*!
    * @brief Construct communication group.
    *
    * The number of children per node is flexible.
    *
    * @param nchild Number of children per node in the group,
    *        i.e., nchild=2 is a binary tree.
    */
   AsyncCommGroup( const int nchild );


   /*!
    * @brief Destructor.
    */
   virtual ~AsyncCommGroup(void);


   //@{
   //! @name Define the communication group

   /*!
    * @brief Setup the tree for the given group of processes.
    * The root process is specified by its index in the group array.
    *
    * The root is specified by dereferencing @c group
    * array with @c root_index.
    */
   void setGroupAndRootIndex( const Array<int> &group,
                              const int root_index );

   /*!
    * @brief Setup the group for the given group.
    * The root process is specified by its rank.
    *
    * The rank of the root is root_rank, which must
    * be one of the ranks given in the group.
    */
   void setGroupAndRootRank( const Array<int> &group,
                             const int root_rank );

   //@}


   /*!
    * @brief Set the MPI tag used for communication within the group.
    *
    * @attention This class is NOT (and cannot be) responsible for
    * ensuring that the MPI communicator and tag are sufficient to
    * select the correct messages.  Please specify appropriate values
    * for the MPI communicator and tag.  Very elusive bugs can occur
    * if incorrect messages are received.
    */
   void setMPITag( const int mpi_tag );


   /*!
    * @brief Set the MPI communicator used for communication within the group.
    *
    * @attention This class is NOT (and cannot be) responsible for
    * ensuring that the MPI communicator and tag are sufficient to
    * select the correct messages.  Please specify appropriate values
    * for the MPI communicator and tag.  Very elusive bugs can occur
    * if incorrect messages are received.  To be safe, it is best
    * to create a new communicator to avoid interference with other
    * communications within SAMRAI.
    */
   void setMPICommunicator( SAMRAI_MPI::comm &mpi_communicator );


   /*!
    * @brief Set whether sends to parents should be blocking.
    *
    * The default is to use blocking send to parent.
    * Because there is just one parent, short messages
    * can be buffered by MPI to improve the performance
    * of blocking sends.  Blocking sends need not be checked
    * for completion.
    */
   void setUseBlockingSendToParent( const bool flag );


   /*!
    * @brief Set whether sends to children should be blocking.
    *
    * The default is to use nonblocking send to children.
    * Nonblocking sends to children are generally appropriate
    * as there are multiple children.
    */
   void setUseBlockingSendToChildren( const bool flag );




   //@{

   /*!
    * @name Communication methods
    */

   /*!
    * @brief Begin a broadcast communication.
    *
    * Root process of broadcast may send less data
    * (smaller size) than receivers of broadcast,
    * in which case the missing data is considered
    * irrelevant by the root.
    *
    * If this method returns false, checkBcast() must
    * be called until it returns true before any change
    * in object state is allowed.
    *
    * @return Whether operation is completed.
    */
   bool beginBcast( int *buffer, int size );

   /*!
    * @brief Check the current broadcast communication
    * and complete the broadcast if all MPI requests
    * are fulfilled.
    *
    *
    * If no communication is in progress, this call
    * does nothing.
    *
    * @return Whether operation is completed.
    */
   bool checkBcast();


   /*!
    * @brief Begin a gather communication.
    *
    * Sending processes of gather may send less data
    * (smaller size) than receivers of broadcast,
    * in which case the missing data is considered
    * irrelevant by the sender.
    *
    * If this method returns false, checkGather() must
    * be called until it returns true before any change
    * in object state is allowed.
    *
    * On non-root processes, buffer should contain
    * the data to be gathered.  On the root process,
    * it should have enough space for all the data
    * from all the processes in the group.
    *
    * @param buffer Data to gather.
    * @param size Size of data contributed by each process.
    *
    * @return Whether operation is completed.
    */
   bool beginGather( int *buffer, int size );

   /*!
    * @brief Check the current gather communication
    * and complete the gather if all MPI requests
    * are fulfilled.
    *
    * @return Whether operation is completed.
    */
   bool checkGather();


   /*!
    * @brief Begin a sum reduce communication.
    *
    * Assume all messages are the same size.
    *
    * If this method returns false, checkSumReduce() must
    * be called until it returns true before any change
    * in object state is allowed.
    *
    * Buffer should contain the data to be gathered.
    *
    * @return Whether operation is completed.
    */
   bool beginSumReduce( int *buffer, int size );

   /*!
    * @brief Check the current sum reduce communication
    * and complete the sum reduce if all MPI requests
    * are fulfilled.
    *
    * @return Whether operation is completed.
    */
   bool checkSumReduce();


   /*!
    * @brief Check the current communication
    * and complete it if all MPI requests
    * are fulfilled.
    */
   bool checkOperation();

   /*!
    * @brief Wait for the current operation to complete.
    */
   void waitOperation();


   /*!
    * @brief Whether the last communication operation has finished.
    */
   bool isDone() const;

   //@}


   int getNumberOfChildren() const;


  void logCurrentState( std::ostream &co ) const;


private:


   /*!
    * @brief Return the requests for use by this object.
    *
    * This object allocates the MPI_Request objects it requires.
    * However, to operate within a communication stage (such as
    * AsyncCommStage), this method may be overiden to provide
    * externally managed requests.
    */
   virtual SAMRAI_MPI::request *getRequestPointer() const;


   /*
    * @brief Assert that user-set MPI parameters are valid.
    */
   void checkMPIParams();



#ifdef HAVE_MPI


   /*!
    * @brief Operation disallowed due to primitive internal memory management.
    */
   AsyncCommGroup(const AsyncCommGroup &r) : d_nchild(0) {};

   /*!
    * @brief Operation disallowed by primitive internal memory management.
    */
   AsyncCommGroup &operator=(const AsyncCommGroup &r) { return *this; };


   //! @brief Operations user would want to do.
   enum BaseOp { undefined,
                 gather,
                 bcast,
                 max_reduce,
                 min_reduce,
                 sum_reduce };
   //! @brief Tasks, executed in order, to complete a base operation.
   enum TaskOp { recv_start,
                 recv_check,
                 send_start,
                 send_check,
                 none };

   struct ChildData {
      //! @brief Rank of child process in the group.
      int rank;
      //! @brief Number of descendants on each child branch.
      int size;
      ChildData() : rank(-1), size(-1) {}
   };


   /*!
    * @brief Begin a generic reduce communication.
    *
    * This method is the workhorse underneath the public reduce methods.
    *
    * If this method returns false, checkOperation() must
    * be called until it returns true before any change
    * in object state is allowed.
    *
    * Buffer should contain the data to be gathered.
    *
    * @return Whether operation is completed.
    */
   bool beginReduce();

   /*!
    * @brief Check the current gather communication.
    *
    * This method is the workhorse underneath the public reduce methods.
    *
    * @return Whether operation is completed.
    */
   bool checkReduce();

   /*!
    * @brief Perform reduction on data that after it has been brought
    * to the local process.
    *
    * The exact reduce operation depends on the base operation.
    */
   void reduceData( int *output, int *data ) const;


   /*!
    * @brief Compute the data that depends on the group definition.
    *
    * Extract and compute parts and characteristics of the tree
    * relevant to the local process.
    */
   void computeDependentData( const Array<int> &group );

   void resetStatus();



   //@{
   /*!
    * @name Mappings between array indices, group positions and process ranks
    */

   /*
    * pos refers the process position in the group (where root position == 0)
    * idx refers the index of the process in the group
    * rank refers the process rank
    */

   /*!
    * @brief Convert the array index to the position.
    */
   int toPosition( int index ) const;
   /*!
    * @brief Convert the position to the array index.
    */
   int toIndex( int position ) const;

   /*!
    * @brief Compute the position of child child_id of a parent (whether
    * or not that child really exists).
    *
    * @param parent_pos Position of the parent in the group.
    * @param ic Index of the child.  (Zero coresponds to the first child.)
    */
   int toChildPosition( int parent_pos, int ic ) const;

   /*!
    * @brief Compute the oldest (lowest position) child position
    * of a given position (whether or not that child really exists).
    *
    * Same as toChildPosition( parent_pos, 0 );
    */
   int toOldest( int parent_pos ) const;
   /*!
    * @brief Compute the youngest (highest position) child position
    * of a given position (whether or not that child really exists).
    *
    * Same as toChildPosition( parent_pos, d_nchild-1 );
    */
   int toYoungest( int parent_pos ) const;

   //@}


   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();


   /*!
    * @brief Number of children per node.
    */
   const int d_nchild;


   /*!
    * @brief Index of the local process in d_group_ranks.
    *
    * The group is defined by an array of process ranks.
    * The index of a process refers to the index in the this array.
    * The "position" of a process in the group represents
    * the position in the group, where the nodes
    * are numbered sequentially, starting at zero for the root.
    *
    * We require that the root has position zero.
    * If the index of the root is zero, then positions are
    * identical to indices.  If not, then the index of the
    * root is swapped with zero to get positions.  These two
    * cases correspond respectively to following trivial and
    * nontrivial maps.
    *
    * @verbatim
    *
    *                 Trivial map           Nontrivial map
    * Parameter     d_root_idx == 0         d_root_idx > 0
    * ---------     ---------------         --------------
    *
    * index of root       0              d_root_idx
    *
    * index of          d_idx            d_idx
    * local process
    *
    * index of            p              p == 0 ? d_root_idx :
    * position p                         p == d_root_idx ? 0 :
    *                                    p
    *
    * position of         i              i == 0 ? d_root_idx :
    * index i                            i == d_root_idx ? 0 :
    *                                    i
    *
    * @endverbatim
    */
   int d_idx;

   /*!
    * @brief Index of the root process in d_group_ranks.
    */
   int d_root_idx;

   int d_group_size;

   /*!
    * @brief Rank of parent process in the group.
    *
    * If negative, there is no parent (this is the root).
    * In send_start tasks, send only to children with valid ranks
    * (not -1).
    */
   int d_parent_rank;


   //! @brief Data on each child branch.
   ChildData *d_child_data;

   /*!
    * @brief Total of all branch sizes.
    */
   int d_branch_size_totl;


   /*!
    * @brief Operation being performed.
    */
   BaseOp d_base_op;

   /*!
    * @brief Next task in a current communication operation.
    *
    * If d_next_task_op is none, there is no current communication
    * operation (the last one is completed).
    */
   TaskOp d_next_task_op;


   /*!
    * @brief External data input and output buffer.
    *
    * This provides the input and output for transfering data.
    * The expected size of the buffer depends on the communication.
    */
   int *d_external_buf;

   /*!
    * @brief Size of d_external_buf.
    */
   int d_external_size;


   /*!
    * @brief Internal buffer.
    *
    * Used for gather and reduce operations but not for
    * broadcast.
    */
   Array<int> d_internal_buf;

   /*!
    * @brief Requests managed internally.
    *
    * This is set to NULL if getRequestPointer() is overriden to
    * provide the requests externally.
    */
   mutable SAMRAI_MPI::request *d_internal_requests;

   int d_mpi_tag;
   SAMRAI_MPI::comm d_mpi_communicator;

   bool d_use_blocking_send_to_children;
   bool d_use_blocking_send_to_parent;




   // Make some temporary variable statuses to avoid repetitious allocations.
   SAMRAI_MPI::status d_mpi_status;

   int d_mpi_err;

   static Pointer<Timer> t_reduce_data;
   static Pointer<Timer> t_wait_all;

#ifdef DEBUG_CHECK_ASSERTIONS
   /*!
    * @brief Array of process ranks in the group.
    *
    * It is possible to code this class without storing all the
    * ranks internally.  However, it is easier to debug if the
    * ranks are available.
    */
   Array<int> d_group_ranks;
#endif


#else


   /*!
    * @brief Operation disallowed due to primitive internal memory management.
    */
   AsyncCommGroup(const AsyncCommGroup &r) : d_nchild(0) {};

   /*!
    * @brief Operation disallowed by primitive internal memory management.
    */
   AsyncCommGroup &operator=(const AsyncCommGroup &r) { return *this; };

   const int d_nchild;

#endif

};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/AsyncCommGroup.I"
#endif

#endif  // included_tbox_AsyncCommGroup
