/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Schedule of communication transactions between processors
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {


/*
 *************************************************************************
 * Access number of send transactions.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
Schedule::getNumLocalTransactions() const
{
   return d_local_set.size();
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::setMPITag(
   const int tag0,
   const int tag1)
{
   TBOX_ASSERT(tag0 >= 0);
   TBOX_ASSERT(tag1 >= 0);
   d_first_tag = tag0;
   d_second_tag = tag1;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::setMPI(
   const SAMRAI_MPI& mpi)
{
   d_mpi = mpi;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::deallocateCommunicationObjects()
{
   delete[] d_coms;
   d_coms = NULL;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::setFirstMessageLength(
   int first_message_length)
{
   TBOX_ASSERT(first_message_length > 0);
   d_first_message_length = first_message_length;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::initializeCallback()
{
   TimerStruct& timers(s_static_timers[s_default_timer_prefix]);
   getAllTimers(s_default_timer_prefix, timers);
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Schedule::finalizeCallback()
{
   s_static_timers.clear();
}

}
}
