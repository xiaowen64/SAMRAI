//
// File:  $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BoxComm.C $
// Package:  SAMRAI hierarchy
// Copyright:  (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:  $LastChangedRevision: 2141 $
// Modified:  $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Utility class for parallel communication of boxes
//

#ifndef included_hier_BoxComm_C
#define included_hier_BoxComm_C

#include "BoxComm.h"



namespace SAMRAI {
   namespace hier {

/*
 * ************************************************************************
 * 
 * Broadcast a Box from one to all processors
 * 
 * ************************************************************************
 */

template<int DIM> void BoxComm<DIM>::bcastBox(Box<DIM> &box, const int root, tbox::SAMRAI_MPI::comm comm)
{
   if (tbox::SAMRAI_MPI::getNodes() > 1) {
      int *buffer = new int[DIM*2];
      for (int i=0; i<DIM; ++i) {
        buffer[i*2] = box.lower(i);
        buffer[i*2+1] = box.upper(i);
      }

      int size = DIM*2;
      tbox::SAMRAI_MPI::comm tmp = tbox::SAMRAI_MPI::getCommunicator();
      tbox::SAMRAI_MPI::setCommunicator(comm);
      tbox::SAMRAI_MPI::bcast(buffer, size, root);
      tbox::SAMRAI_MPI::setCommunicator(tmp);

      for (int j=0; j<DIM; ++j) {
         box.lower(j) = buffer[j*2];
         box.upper(j) = buffer[j*2+1];
      }
      delete [] buffer;
   }
}

/*
 * ************************************************************************
 * 
 * Broadcast a BoxList from one to all processors
 * 
 * ************************************************************************
 */

template<int DIM> void BoxComm<DIM>::bcastBoxList(BoxList<DIM> &box_list, const int root)
{
#ifdef HAVE_MPI
   BoxArray<DIM> tmp(box_list);
   bcastBoxArray(tmp, root);
   box_list = tmp;
#endif
}

/*
 * ************************************************************************
 * 
 * Broadcast a BoxArray from one to all processors
 * 
 * ************************************************************************
 */

template<int DIM> void BoxComm<DIM>::bcastBoxArray(BoxArray<DIM> &box_array, const int root)
{
   if (tbox::SAMRAI_MPI::getNodes() > 1) {
      /*
       * root broadcasts the number of boxes in the list
       */
      int len = 0;
      if (tbox::SAMRAI_MPI::getRank() == root) {
         len = box_array.getNumberOfBoxes();
      }

      len = tbox::SAMRAI_MPI::bcast(len, root);

      box_array.resizeBoxArray(len);

      /*
       * If the box_list is empty, don't need to do anything!
       */
      if (len != 0) {
   
         /*
          * allocate buffer to hold box coordinates
          */
         int *buffer = new int[len*DIM*2];
   
         /*
          * root fills the buffer for sending the boxes' coordinates
          * to all other processors.
          */
         if (tbox::SAMRAI_MPI::getRank() == root) {

            for (int b=0; b<len; ++b) {
               int offset  = b*DIM*2;
               for (int i=0; i<DIM; ++i) {
                  buffer[offset++] = box_array[b].lower(i);
                  buffer[offset++] = box_array[b].upper(i);
               }
            }
         }
   
         /*
          * broadcast box coordinates to all processors
          */
         int buf_len = len*DIM*2;
         buf_len = tbox::SAMRAI_MPI::bcast(buf_len, root);
         tbox::SAMRAI_MPI::bcast(buffer, buf_len, root);
   
         /*
          * all processors except root (whose data we promised to treat
          * as const) construct box_array.
          */
         if (tbox::SAMRAI_MPI::getRank() != root) {
            box_array.resizeBoxArray(len);
            int offset  = 0;
            for (int b=0; b<len; ++b) {
               for (int j=0; j<DIM; ++j) {
                  box_array[b].lower(j) = buffer[offset++];
                  box_array[b].upper(j) = buffer[offset++];
               }
            }
         }
   
         delete [] buffer;
      }
   }
}

/*
 * ************************************************************************
 * 
 * Send a Box from this processor to another processor.
 * 
 * ************************************************************************
 */
template<int DIM> void BoxComm<DIM>::sendBox(const Box<DIM> &box, const int rcvr_id)
{
   int *buffer = new int[DIM*2];
   for (int i=0; i<DIM; ++i) {
      buffer[i*2] = box.lower(i);
      buffer[i*2+1] = box.upper(i);
   }
   tbox::SAMRAI_MPI::send(buffer, DIM*2, rcvr_id, false);
   delete [] buffer;
}

/*
 * ************************************************************************
 * 
 * Send a Box from this processor to several processors.
 * 
 * ************************************************************************
 */

template<int DIM> void BoxComm<DIM>::sendBox(const Box<DIM> &box, 
                            const tbox::Array<int> &rcvr_id)
{
#ifdef HAVE_MPI
   int ct = rcvr_id.getSize();
   MPI_Request *request = new MPI_Request[ct];
   MPI_Status *status = new MPI_Status[ct];
   
   int *buffer = new int[DIM*2];
   for (int i=0; i<DIM; ++i) {
      buffer[i*2] = box.lower(i);
      buffer[i*2+1] = box.upper(i);
   }
   
   //start non-blocking sends to each descendant
   for (int j=0; j<ct; ++j) {
      int dest = rcvr_id[j];
      MPI_Isend(buffer, DIM*2, MPI_INT, dest,
    0, MPI_COMM_WORLD, request+j);
   }
   
   //wait for all sends to finish
   MPI_Waitall(ct, request, status);
   delete [] request;
   delete [] status;

   delete [] buffer;
#endif
}


/*
 * ************************************************************************
 * 
 * Receive a Box from another processor; must be
 * used in conjunction with the matching call, sendBox();
 * 
 * ************************************************************************
 */
template<int DIM> void BoxComm<DIM>::recvBox(Box<DIM> &box, const int sender_id)
{
   int size = DIM*2;
   int *buffer = new int[size];
      
   tbox::SAMRAI_MPI::recv(buffer, size, sender_id, false);
   for (int j=0; j<DIM; ++j) {
      box.lower(j) = buffer[j*2];
      box.upper(j) = buffer[j*2+1];
   }
   delete [] buffer;
}

/*
 * ************************************************************************
 * 
 * all-to-all exchange of box arrays and associated weights
 * 
 * ************************************************************************
 */
template<int DIM> void BoxComm<DIM>::exchangeBoxArraysAndWeightArrays(
   const BoxArray<DIM> &box_array_in,
   const tbox::Array<double> &weights_in,
   BoxArray<DIM> &box_array_out,
   tbox::Array<double> &weights_out)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(box_array_in.getNumberOfBoxes() == weights_in.getSize());
#endif

   /*
    * allocate send and receive buffers, and set array sizes 
    * for the output arrays.
    */
   int size_in = box_array_in.getNumberOfBoxes();
   int size_out = tbox::SAMRAI_MPI::sumReduction(size_in);

#ifdef DEBUG_CHECK_ASSERTIONS
   if (size_out <= 0) {
      TBOX_ERROR("BoxComm<DIM>::exchangeBoxArraysAndWeightArrays() error"
                 << "\n All box arrays have zero length!" << std::endl);
   }
#endif

   int buf_size_in  = size_in*DIM*2;
   int buf_size_out = size_out*DIM*2;
 
   box_array_out.resizeBoxArray(size_out);
   weights_out.resizeArray(size_out);
 
   tbox::Array<int> buf_in(buf_size_in);
   tbox::Array<int> buf_out(buf_size_out);
 
   int* buf_in_ptr = (int*)NULL;
   int* buf_out_ptr = (int*)NULL;
   const double* wgts_in_ptr = (const double*)NULL;
   double* wgts_out_ptr = (double*)NULL;

   if (size_in > 0) {
      buf_in_ptr = buf_in.getPointer(); 
      wgts_in_ptr  = weights_in.getPointer();
   }
   if (size_out > 0) {
      wgts_out_ptr = weights_out.getPointer();
      buf_out_ptr = buf_out.getPointer(); 
   }
 
   /*
    * populate the buffers with data for sending
    */
   int offset = 0;
   for (int x = 0; x < size_in; ++x) {
     for (int i = 0; i < DIM; ++i) {
	buf_in_ptr[offset++] = box_array_in[x].lower(i);
	buf_in_ptr[offset++] = box_array_in[x].upper(i);
     }
   }
 
   /*
    * exchange the data 
    */
   tbox::SAMRAI_MPI::allGather(buf_in_ptr, buf_size_in, buf_out_ptr, buf_size_out);
   tbox::SAMRAI_MPI::allGather(wgts_in_ptr, size_in, wgts_out_ptr, size_out);
 
   /*
    * assemble the output array of boxes
    */
   offset  = 0;
   for (int b = 0; b < size_out; ++b) {
     for (int j = 0; j < DIM; ++j) {
       box_array_out[b].lower(j) = buf_out_ptr[offset++];
       box_array_out[b].upper(j) = buf_out_ptr[offset++];
     }
   }

}
 

}
}

#endif
