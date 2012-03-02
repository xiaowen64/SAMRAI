/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   All-to-one and one-to-all communication using a tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::setMPITag(
   const int mpi_tag)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI tag is not allowed\n"
         << "during pending communications");
   }
   d_mpi_tag = mpi_tag;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::setUseMPICollectiveForFullGroups(
   bool use_mpi_collective)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI collective option is not allowed\n"
         << "during pending communications");
   }
   d_use_mpi_collective_for_full_groups = use_mpi_collective;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::setUseBlockingSendToChildren(
   const bool flag)
{
   d_use_blocking_send_to_children = flag;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::setUseBlockingSendToParent(
   const bool flag)
{
   d_use_blocking_send_to_parent = flag;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::getNumberOfChildren() const
{
   return static_cast<int>(d_nchild);
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::toPosition(
   int index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (index < 0 || index >= d_group_size) {
      TBOX_ERROR(
         "Invalid index " << index << "\n"
                          << "should be in [0," << d_group_size - 1
                          << "].");
   }
#endif
   const int position = index == 0 ? d_root_idx :
      index == d_root_idx ? 0 :
      index;
   return position;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::toIndex(
   int position) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (position < 0 || position >= d_group_size) {
      TBOX_ERROR(
         "Invalid parent position " << position
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   const int index = position == 0 ? d_root_idx :
      position == d_root_idx ? 0 :
      position;
   return index;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::toChildPosition(
   int parent_pos,
   int child) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= d_group_size) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return parent_pos * static_cast<int>(d_nchild) + 1 + child;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::toOldest(
   int parent_pos) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= d_group_size) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return parent_pos * static_cast<int>(d_nchild) + 1;
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
AsyncCommGroup::toYoungest(
   int parent_pos) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= d_group_size) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return (parent_pos + 1) * static_cast<int>(d_nchild);
}

/*
 ****************************************************************
 ****************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::resetStatus()
{
   /*
    * Do not set d_mpi_status.count because it is not guaranteed
    * to exist in every implemenation.
    */
   d_mpi_status.MPI_TAG =
      d_mpi_status.MPI_SOURCE =
         d_mpi_status.MPI_ERROR = -1;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
AsyncCommGroup::bcastByMpiCollective()
{
   d_mpi.Bcast(d_external_buf,
      d_external_size,
      MPI_INT,
      d_root_rank);
   d_next_task_op = none;
   return true;
}

/*
 **********************************************************************
 * The fact that the current reduction is a sum reduction is recorded
 * in the d_base_op parameter.  We simply use the generic method for
 * checking reduction.
 **********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
AsyncCommGroup::checkSumReduce()
{
   return checkReduce();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::checkMPIParams()
{
   if (d_mpi_tag < 0) {
      TBOX_ERROR("AsyncCommGroup: Invalid MPI tag value "
         << d_mpi_tag << "\nUse setMPITag() to set it.");
   }
}

/*
 ***************************************************************************
 * Allocate static timers.
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::initializeCallback()
{
   t_reduce_data = TimerManager::getManager()->
      getTimer("tbox::AsyncCommGroup::reduceData()");
   t_wait_all = TimerManager::getManager()->
      getTimer("tbox::AsyncCommGroup::mpi_wait_all");
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
AsyncCommGroup::finalizeCallback()
{
   t_reduce_data.reset();
   t_wait_all.reset();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
