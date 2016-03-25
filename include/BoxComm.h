//
// File:  $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxComm.h $
// Package:  SAMRAI hierarchy
// Copyright:  (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:  $LastChangedRevision: 1704 $
// Modified:  $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:  Utility class for parallel communication of boxes
//

#ifndef included_hier_BoxComm
#define included_hier_BoxComm

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_SAMRAI_MPI
#include "tbox/SAMRAI_MPI.h"
#endif


namespace SAMRAI {
   namespace hier {


/**
 * Class BoxComm<DIM> is a utility class that provides support for
 * for broadcast and point-to-point communication of Boxes, BoxLists,
 * and BoxArrays.  All calls in this class are static.  Point-to-point
 * communication is blocking.  All functions fall through when SAMRAI
 * is compiled without MPI.
 * 
 * Note that this class is a utility class to group function calls in one
 * name space (all calls are to static functions).  Thus, you should never
 * attempt to instantiate a class of type BoxComm<DIM>; simply call the
 * functions as static functions using the BoxComm<DIM>::function(...)
 * syntax.
 * 
 * @see hier::Box
 * @see hier::BoxArray
 * @see hier::BoxList
 */

template<int DIM> struct BoxComm
{
   /**
    * Broadcast a Box from specified root process to all other processors
    * (root's box is treated as const).
    */
   static void bcastBox(Box<DIM> &box, 
                        const int root=0, 
                        tbox::SAMRAI_MPI::comm = tbox::SAMRAI_MPI::commWorld);

   /**
    * Broadcast a BoxList from specified root process to all other processors.
    * Processors other than root do NOT need to know the number of boxes in
    * the list prior to calling, since box_list will be re-sized to
    * the size of root's list (root's list is treated as const).
    */
   static void bcastBoxList(BoxList<DIM> &box_list, const int root=0);

   /**
    * Broadcast a BoxArray from specified root process to all other processes.
    * Processors other than root do NOT need to know the number of boxes in
    * the array prior to calling, since box_array will be re-sized to
    * the size of root's array (root's array is treated as const).
    */
   static void bcastBoxArray(BoxArray<DIM> &box_array, const int root=0);

   /**
    * Send a Box from this processor to another processor.
    * Blocking communications are used.  This call must be matched
    * by a call to recvBox, or your code will hang.
    */
   static void sendBox(const Box<DIM> &box, const int rcvr_id);

   /**
    * Send a Box from this processor to several processors.
    */
   static void sendBox(const Box<DIM> &box, const tbox::Array<int> &rcvr_id);

   /**
    * Receive a Box from the designated processor (sender_id).
    * Blocking communications are used.  This call must be matched
    * by a call to sendBox, or your code will hang.
    */
   static void recvBox(Box<DIM> &box, const int sender_id);

   /**
    * All-to-all communication of box arrays and associated weights.
    * On invocation, each processor has a (possibly empty) array of 
    * 'owned' boxes, and each box has a weight.  On return, 
    * each processor has a array that contains all boxes owned 
    * by all processors, and their associated weights.
    * If all processors input arrays have zero length, an error
    * is thrown.
    */
   static void exchangeBoxArraysAndWeightArrays(
      const BoxArray<DIM> &box_array_in,
      const tbox::Array<double> &weights_in,
      BoxArray<DIM> &box_array_out,
      tbox::Array<double> &weights_out);
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxComm.C"
#endif

#endif

