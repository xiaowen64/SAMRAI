/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Staged peer-to-peer communication.
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
template<class TYPE>
SAMRAI_INLINE_KEYWORD
void AsyncCommPeer<TYPE>::setMPITag(
   const int tag0,
   const int tag1)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI tag is not allowed\n"
         << "during pending communications");
   }
   d_tag0 = tag0;
   d_tag1 = tag1;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
void AsyncCommPeer<TYPE>::setMPI(
   const SAMRAI_MPI& mpi)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI object is not allowed\n"
         << "during pending communications");
   }
   d_mpi = mpi;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool AsyncCommPeer<TYPE>::isDone() const
{
   return d_next_task_op == none;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
void AsyncCommPeer<TYPE>::setPeerRank(
   int peer_rank)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the peer is not allowed\n"
         << "during pending communications");
   }
   d_peer_rank = peer_rank;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
int AsyncCommPeer<TYPE>::getPeerRank() const
{
   return d_peer_rank;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
void AsyncCommPeer<TYPE>::limitFirstDataLength(
   size_t limit)
{
   d_max_first_data_len = limit;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
void AsyncCommPeer<TYPE>::resetStatus(
   SAMRAI_MPI::Status& mpi_status)
{
   /*
    * Do not set stat.count because it is not guaranteed
    * to exist in every implemenation.
    */
   mpi_status.MPI_TAG =
      mpi_status.MPI_SOURCE =
         mpi_status.MPI_ERROR = -1;
}

/*
 ****************************************************************
 ****************************************************************
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
size_t AsyncCommPeer<TYPE>::getNumberOfFlexData(
   size_t number_of_type_data) const
{
   size_t number_of_flexdata = number_of_type_data * sizeof(TYPE);
   number_of_flexdata = number_of_flexdata / sizeof(FlexData)
      + (number_of_flexdata % sizeof(FlexData) > 0);
   return number_of_flexdata;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool AsyncCommPeer<TYPE>::initialize(
   void)
{
   static StartupShutdownManager::Handler handler(
      0,
      0,
      0,
      AsyncCommPeer::finalizeCallback,
      tbox::StartupShutdownManager::priorityTimers);

   t_waitall_timer = TimerManager::getManager()->
      getTimer("tbox::AsyncCommPeer::wait_all()");
   t_send_timer = TimerManager::getManager()->
      getTimer("tbox::AsyncCommPeer::MPI_ISend");
   t_recv_timer = TimerManager::getManager()->
      getTimer("tbox::AsyncCommPeer::MPI_Irecv");

   return true;
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
