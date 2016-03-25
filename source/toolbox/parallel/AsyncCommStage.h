/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/parallel/AsyncCommStage.h $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2132 $
 * Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
 * Description: Support for coordinating multiple asynchronous communications
 */

#ifndef included_tbox_AsyncCommStage
#define included_tbox_AsyncCommStage

#include "SAMRAI_config.h"

#include "tbox/Array.h"

#include "tbox/RelaunchableJob.h"

#include "tbox/AsyncCommGroup.h"

#include "tbox/SAMRAI_MPI.h"

#include "tbox/Timer.h"

namespace SAMRAI {
   namespace tbox {

/*!
 * @brief Stage multiple asynchronous group communications
 * so that the collective can advance asynchronously
 * (as individual underlying MPI requests are completed).
 *
 * Use this class when you
 * - Have multiple group communications, each performed
 *   by a different AsyncCommGroup object.
 * - Want the multiple communications to proceed
 *   independently and asynchronously.
 *
 * This class allocates a set of AsyncCommGroup objects and
 * manages the space for their SAMRAI_MPI::request's.  The requests are
 * staged such that a single MPI_Waitany or MPI_Waitsome call
 * applies to all communication groups allocated by the stage.
 * Thus the communications performed by the AsyncCommGroup
 * objects can complete in the order allowed by the MPI messages.
 * The exact order is NOT deterministic!
 *
 * To advance the communication operation of any of the allocated
 * AsyncCommGroup objects, use advanceAny() or advanceSome().
 * In general, advanceSome() has better performance than advanceAny()
 * because it gets around the "starvation" problem.  See the MPI
 * documentation for a discussion of starvation.
 *
 * This class supports communication and uses MPI for
 * message passing.  If MPI is disabled, the job of this
 * class disappears and the class is effectively empty,
 * except for allocating and deallocating AsyncCommGroup
 * objects.  
 * The public interfaces still remain so the class
 * can compile, but the implementations are trivial.
 */
class AsyncCommStage
{

public:


   /*!
    * @brief Construct a stage that may begin allocating and
    * managing groups.
    */
   AsyncCommStage();


   /*!
    * @brief Deallocate groups remaining in the stage and all
    * internal data used to manage groups.
    */
   virtual ~AsyncCommStage(void);


   /*!
    * @brief Allocate a group in this stage.
    *
    * You may specify a handle associated with the allocated group
    * and pointing back to an object.
    * These pointers are returned when their group communications
    * are completed by the stage methods advanceAny() or advanceSome().
    * They help identify the user object waiting for the communication
    * to complete.
    *
    * To use the handle, create a user class inheritting from
    * RelaunchableJob so that any handle can be dynamically
    * casted back to the user class.
    *
    * Usage outline:
    * -# Derive class UserClass from RelaunchableJob.
    *    Implement the pure virtual interfaces therein.
    * -# Associate a UserClass object with a AsyncCommGroup
    *    object by using the address of the UserClass object to
    *    allocate the group.  See AsyncCommStage::allocateCommGroup().
    * -# Use the allocated group to initiate a group communication.
    * -# Repeat the above for any number of user objects.
    * -# Determine the specific user object associated with the
    *    communication completed by AsyncCommStage::advanceAny()
    *    and AsyncCommStage::advanceSome() by dynamically
    *    casting the output RelaunchableJob pointer to 
    *    pointers to UserClass objects.
    *
    * Usage example:
    * @verbatim
    * class UserClass : RelaunchableJob { ... };
    * UserClass user_object[10];
    * AsyncCommStage stage;
    * for ( i=0; i<10; ++i ) {
    *    AsyncCommGroup *group =
    *      stage.allocateCommGroup( 2, &user_object[i] );
    *    useGroupToInitiateCommunication(group);
    * }
    * RelaunchableJob *handle = stage.advanceAny();
    * UserClass *completed_user_object = dynamic_cast<UserClass*>(handle);
    * @endverbatim
    *
    * Allocating a group with a NULL handle works, but
    * you lose the ability to map a specific communication
    * back to a user object.  See advanceAny() and advanceSome().
    *
    * An allocated group should be deallocated using the delete operator
    * when no longer needed.  All allocated groups will be deleted
    * when the stage goes out of scope.  Since it is an error to
    * delete a group that has pending message requests, it is also
    * an error to delete a stage that has groups with pending
    * message requests.
    *
    * @param nchild Number of child branches per tree node
    *               (see AsyncCommGroup).
    * @param handle Optional pointer back to a user object associated
    *               with the group.  See class documentation.
    */
   AsyncCommGroup *allocateCommGroup(
      int nchild,
      RelaunchableJob *handle=NULL );


   /*!
    * @brief Advance to completion one group (any group) that is
    * currently waiting for communication to complete.
    *
    * This method uses MPI_Waitany, which may be prone to starvation.
    * It is better to use advanceSome(), which uses MPI_Waitsome,
    * which avoids starvation.
    *
    * @return the RelaunchableJob pointer used to
    * construct the group that is completed.  NULL if no
    * group in the stage had incomplete communication.
    */
   RelaunchableJob *advanceAny();


   /*!
    * @brief Advance to completion one or more groups (any groups)
    * that are currently waiting for communication to complete.
    *
    * @param completed Array space to put handleegies associated
    *                  with the completed communications.
    *
    * @return number of groups completed or 0 if none have outstanding
    * communication.  The strategy pointers corresponding to the
    * completed groups are set in the @c completed array.
    */
   int advanceSome( Array<RelaunchableJob*> &completed );



   /*!
    * @brief Return whether the stage is done with all requested
    * communications.
    */
   bool isDone() const;

   /*!
    * @brief Return the number of allocated groups that have pending
    * communication.
    */
   bool numberOfOutstandingGroups() const;

   /*!
    * @brief Return the number of pending SAMRAI_MPI::request objects.
    */
   bool numberOfOutstandingRequests() const;



private:


   /*!
    * @brief Augmented AsyncCommGroup, to operate on a stage.
    *
    * This class adds data to AsyncCommGroup so it can be
    * managed on a stage and access requests managed on the
    * same stage.
    */
   struct StagedGroup : public AsyncCommGroup {
     StagedGroup( const int nchild,
                  AsyncCommStage *stage,
                  const int index,
                  RelaunchableJob *handle );
     ~StagedGroup();
     /*!
      * @brief The stage that generated this object
      * and can provide SAMRAI_MPI::request objects for it.
      */
     AsyncCommStage *d_stage;

     //! @brief Group's index within the stage, assigned by the stage.
     const int d_index;

     /*!
      * @brief Handle associated with group, used by the stage.
      *
      * Currently, none of the implementations required by RelaunchableJob
      * is used by AsyncCommStage, so it is possible to have trivial
      * implementations of the RelaunchableJob base class.
      * However, this may change in the future, or you may be
      * using other classes that accesses the RelaunchableJob
      * implementations after they are returned from AsyncCommStage.
      */
     RelaunchableJob *d_handle;

     //!@ Override the virtual method to return requests on the stage.
     SAMRAI_MPI::request *getRequestPointer() const;
   };


   /*
    * Member class StagedGroup is a friend so it can access
    * private methods destageGroup and lookupRequestPointer.
    * This avoids making those private methods public.
    * We are not breaking encapsulation because StagedGroup
    * is a private class.
    */
   friend class StagedGroup;


   /*!
    * @brief Trivial typedef used by SAMRAI template scripts to
    * generate code to instantiate Array<StagedGroup*> (not
    * used anywhere else).
    */
   typedef StagedGroup* StagedGroupPtr;


   /*!
    * @brief Clear references to a group allocated on this stage.
    *
    * To ONLY be be used just before freeing an allocated group.
    * This clears the references so that the group can be
    * deallocated without leaving behind broken references.
    */
   void destageGroup( StagedGroup *group );


#ifdef HAVE_MPI


   /*!
    * @brief Lookup and return the request pointer from the stage for
    * the given group.
    *
    * The given group MUST have been allocated by the stage.
    * The number of requests that the group may use is the
    * same as the number of children the group has on each branch.
    *
    * The pointer is NOT guaranteed to be the same for the life
    * of the group, as a stage may rearange the array of
    * SAMRAI_MPI::request objects.  However, the pointer is valid
    * until the next call to allocateCommGroup().
    *
    * This is private because only friend class StagedGroup should
    * use it.
    */
   SAMRAI_MPI::request *lookupRequestPointer( const int igroup )const;


   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();



   /*!
    * @brief Number of groups (including deallocated groups) that are
    * still occupying space in the stage.
    */
   int d_n_group;

   /*!
    * @brief Number of request slots, including unused ones
    * that belonged to deallocated groups.
    */
   int d_n_req;

   //! @brief Groups managed in this stage.
   Array<StagedGroup*> d_group;

   /*!
    * @brief SAMRAI_MPI::request objects used by the groups.
    *
    * This is mutable because the const method getRequestPointer()
    * needs to return a non-const SAMRAI_MPI::request pointer.
    * The pointer must be non-const for use in MPI calls.
    * getRequestPointer() should be const because no group
    * should require a non-const stage just to get the
    * request allocated for it.
    */
   mutable Array<SAMRAI_MPI::request> d_req;

   //!@brief Map from request index to group index.
   Array<int> d_req_to_group;

   //!@brief Map from group index to (the group's first) request index.
   Array<int> d_group_to_req;

   static Pointer<Timer> t_wait_any;
   static Pointer<Timer> t_wait_some;


#else

   /*
    * Without MPI, this class manages the communication groups,
    * but not their SAMRAI_MPI::request.
    */
   int d_n_group;
   Array<StagedGroup*> d_group;

#endif


};

}
}

#endif  // included_tbox_AsyncCommStage
