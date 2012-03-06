/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Support for coordinating multiple asynchronous communications
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommStage::Member::setHandler(
   Handler* handler)
{
   d_handler = handler;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
AsyncCommStage::Handler*
AsyncCommStage::Member::getHandler() const
{
   return d_handler;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
size_t
AsyncCommStage::Member::numberOfRequests() const
{
   return d_nreq;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::Request*
AsyncCommStage::Member::getRequestPointer() const
{
   if (d_stage == NULL) {
      TBOX_ERROR("AssyncCommStage::Member::getRequestPointer():\n"
         << "Empty stage encountered!\n"
         << "This probably means that the stage that allocated\n"
         << "your AsyncCommGroup has been deallocated.\n"
         << "(and your AssyncCommGroup deallocated too!).\n"
         << "It is an error to deallocate a stage and still\n"
         << "use the Member it allocated.\n");

   }
   return d_stage->lookupRequestPointer(d_index_on_stage);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::Status*
AsyncCommStage::Member::getStatusPointer() const
{
   if (d_stage == NULL) {
      TBOX_ERROR("AssyncCommStage::Member::getStatusPointer():\n"
         << "Empty stage encountered!\n"
         << "This probably means that the stage that allocated\n"
         << "your AsyncCommGroup has been deallocated.\n"
         << "(and your AssyncCommGroup deallocated too!).\n"
         << "It is an error to deallocate a stage and still\n"
         << "use the Member it allocated.\n");

   }
   return d_stage->lookupStatusPointer(d_index_on_stage);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void AsyncCommStage::clearCompletionQueue()
{
   d_completed_members.clear();
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void AsyncCommStage::setCommunicationWaitTimer(
   const boost::shared_ptr<Timer>& communication_timer)
{
   d_communication_timer = communication_timer;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
size_t
AsyncCommStage::numberOfMembers() const
{
   return d_member_count;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
size_t
AsyncCommStage::numberOfRequests(
   size_t index_on_stage) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(index_on_stage < d_members.size());
   TBOX_ASSERT(d_members[index_on_stage] != NULL);
#endif

   const int init_req = static_cast<int>(d_member_to_req[index_on_stage]);
   const int term_req = static_cast<int>(d_member_to_req[index_on_stage + 1]);
   return term_req - init_req;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
AsyncCommStage::hasPendingRequests() const
{
   bool hasPending = false;
   for (size_t ireq = 0; ireq < d_member_to_req[d_members.size()]; ++ireq) {
      if (d_req[ireq] != MPI_REQUEST_NULL) {
         hasPending = true;
         break;
      }
   }
   return hasPending;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
size_t
AsyncCommStage::numberOfPendingRequests() const
{
   size_t npending = 0;
   for (size_t ireq = 0; ireq < d_member_to_req[d_members.size()]; ++ireq) {
      if (d_req[ireq] != MPI_REQUEST_NULL) {
         ++npending;
      }
   }
   return npending;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
size_t
AsyncCommStage::numberOfPendingMembers() const
{
   size_t nmember = 0;
   for (size_t imember = 0; imember < d_members.size(); ++imember) {
      if (d_members[imember] != NULL &&
          d_members[imember]->hasPendingRequests()) {
         ++nmember;
      }
   }
   return nmember;
}

/*
 ****************************************************************
 * Return the request pointer for a communication Member
 * allocated by this object.
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::Request*
AsyncCommStage::lookupRequestPointer(
   const size_t imember) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(imember < d_members.size());
   TBOX_ASSERT(d_members[imember] != NULL);
#endif
   return &d_req[d_member_to_req[imember]];
}

/*
 ****************************************************************
 * Return the status pointer for a communication Member
 * allocated by this object.
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::Status*
AsyncCommStage::lookupStatusPointer(
   const size_t imember) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(imember < d_members.size());
   TBOX_ASSERT(d_members[imember] != NULL);
#endif
   return &d_stat[d_member_to_req[imember]];
}

}
}
