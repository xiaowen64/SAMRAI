/*
 * File:        AsyncCommGroup.I
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2003 The Regents of the University of California
 * Revision:    $Revision: 453 $
 * Modified:    $Date: 2005-06-16 10:19:28 -0700 (Thu, 16 Jun 2005) $
 * Description: All-to-one and one-to-all communication using a tree.
 */

#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif



namespace SAMRAI {
   namespace tbox {



#ifdef HAVE_MPI



inline
void AsyncCommGroup::reduceData( int *output,
                                     int *data ) const
{
   int i;
   switch ( d_base_op ) {
   case max_reduce:
      for ( i=0; i<d_external_size; ++i ) {
         if ( output[i] < data[i] ) output[i] = data[i];
      }
      break;
   case min_reduce:
      for ( i=0; i<d_external_size; ++i ) {
         if ( output[i] > data[i] ) output[i] = data[i];
      }
      break;
   case sum_reduce:
      for ( i=0; i<d_external_size; ++i ) {
         output[i] = output[i] + data[i];
      }
      break;
   default:
      TBOX_ERROR("Library error: d_base_op is somehow corrupted. ");
   }
   return;
}



inline
void AsyncCommGroup::setMPITag( const int mpi_tag )
{
   if ( ! isDone() ) {
      TBOX_ERROR("Resetting the MPI tag is not allowed\n"
                 <<"during pending communications");
   }
   d_mpi_tag = mpi_tag;
   return;
}



inline
void AsyncCommGroup::setMPICommunicator(
   MPI_Comm &mpi_communicator )
{
   if ( ! isDone() ) {
      TBOX_ERROR("Resetting the MPI communicator is not allowed\n"
                 <<"during pending communications");
   }
   d_mpi_communicator = mpi_communicator;
   return;
}



inline
void AsyncCommGroup::setUseBlockingSendToChildren(
   const bool flag )
{
   d_use_blocking_send_to_children = flag;
   return;
}



inline
void AsyncCommGroup::setUseBlockingSendToParent(
   const bool flag )
{
   d_use_blocking_send_to_parent = flag;
   return;
}


inline
bool AsyncCommGroup::isDone() const
{
   return d_next_task_op == none;
}


inline
int AsyncCommGroup::getNumberOfChildren() const
{
   return d_nchild;
}



inline
int AsyncCommGroup::toPosition( int index ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( index < 0 || index >= d_group_size ) {
      TBOX_ERROR("Invalid index " << index << "\n"
                 <<"should be in [0," << d_group_size-1 << "].");
   }
#endif
   const int position = index == 0 ? d_root_idx :
                        index == d_root_idx ? 0 :
                        index;
   return position;
}


inline
int AsyncCommGroup::toIndex( int position ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( position < 0 || position >= d_group_size ) {
      TBOX_ERROR("Invalid parent position " << position
                 <<" should be in [" << 0 << ',' << d_group_size-1 << "].");
   }
#endif
   const int index = position == 0 ? d_root_idx :
                     position == d_root_idx ? 0 :
                     position;
   return index;
}


inline
int AsyncCommGroup::toChildPosition( int parent_pos, int child ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( parent_pos < 0 || parent_pos >= d_group_size ) {
      TBOX_ERROR("Invalid parent position " << parent_pos
                 <<" should be in [" << 0 << ',' << d_group_size-1 << "].");
   }
#endif
   return parent_pos*d_nchild + 1 + child;
}


inline
int AsyncCommGroup::toOldest( int parent_pos ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( parent_pos < 0 || parent_pos >= d_group_size ) {
      TBOX_ERROR("Invalid parent position " << parent_pos
                 <<" should be in [" << 0 << ',' << d_group_size-1 << "].");
   }
#endif
   return parent_pos*d_nchild + 1;
}


inline
int AsyncCommGroup::toYoungest( int parent_pos ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( parent_pos < 0 || parent_pos >= d_group_size ) {
      TBOX_ERROR("Invalid parent position " << parent_pos
                 <<" should be in [" << 0 << ',' << d_group_size-1 << "].");
   }
#endif
   return (parent_pos+1)*d_nchild;
}


inline
void AsyncCommGroup::resetStatus()
{
   /*
    * Do not set d_mpi_status.count because it is not guaranteed
    * to exist in every implemenation.
    */
   d_mpi_status.MPI_TAG =
   d_mpi_status.MPI_SOURCE =
   d_mpi_status.MPI_ERROR = -1;
   return;
}




#else



inline
void AsyncCommGroup::setMPITag( const int mpi_tag )
{
   return;
}



inline
void AsyncCommGroup::setMPICommunicator(
   MPI::comm &mpi_communicator )
{
   return;
}



inline
void AsyncCommGroup::setUseBlockingSendToChildren(
   const bool flag )
{
   return;
}



inline
void AsyncCommGroup::setUseBlockingSendToParent(
   const bool flag )
{
   return;
}


inline
bool AsyncCommGroup::isDone() const
{
   return true;
}


inline
int AsyncCommGroup::getNumberOfChildren() const
{
   return d_nchild;
}




#endif




#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
