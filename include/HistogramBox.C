//
// File:        HistogramBox.C
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 602 $
// Modified:    $Date: 2005-09-06 11:51:31 -0700 (Tue, 06 Sep 2005) $
// Description: Histogram class for computing tagged cell signatures.
//

#ifndef included_mesh_HistogramBox_C
#define included_mesh_HistogramBox_C

#include "HistogramBox.h"

#include "BinaryTree.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "CellData.h"
#include "CellIterator.h"
#include "tbox/MPI.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "HistogramBox.I"
#endif
namespace SAMRAI {
    namespace mesh {

/*
*************************************************************************
*                                                                       *
* Initialization for static data members.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer<tbox::Timer> HistogramBox<DIM>::t_compute;
template<int DIM>
tbox::Pointer<tbox::Timer> HistogramBox<DIM>::t_commwait;

/*
****************************************************************************
*                                                                          *
* Constructors and destructor for HistogramBox<DIM>.                      *
*                                                                          *
****************************************************************************
*/

template<int DIM>  HistogramBox<DIM>::HistogramBox() 
: d_communication_mode(ORIGINAL_MODE)
{
   if ( t_compute.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( t_commwait.isNull() );
#endif
      // Note: We use BergerRigoutsos timers!
      t_compute = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::compute");
      t_commwait = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::commwait");
      tbox::ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }
}

 
template<int DIM>  HistogramBox<DIM>::HistogramBox(const hier::Box<DIM>& box)
: d_communication_mode(ORIGINAL_MODE)
{
   setBox(box);
   resetHistogram();
   if ( t_compute.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( t_commwait.isNull() );
#endif
      // Note: We use BergerRigoutsos timers!
      t_compute = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::compute");
      t_commwait = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::commwait");
      tbox::ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }
}
 
template<int DIM>  HistogramBox<DIM>::~HistogramBox()
{
}

/*
****************************************************************************
*                                                                          *
* Accessory routines for HistogramBox<DIM>.                               *
*                                                                          *
****************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::setBox(const hier::Box<DIM>& box)
{
   d_box = box;
   for (int id = 0; id < DIM; id++) {
     d_histogram[id].resizeArray(d_box.numberCells(id));
   }
}
 
template<int DIM> void HistogramBox<DIM>::resetHistogram()
{
   for (int id = 0; id < DIM; id++) {
      int hi = d_box.numberCells(id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_histogram[id].getSize() == hi);
#endif
      for (int ic = 0; ic < hi; ic++) {
         d_histogram[id][ic] = 0;
      }
   }
}

/*
****************************************************************************
*                                                                          *
* Compute smallest box bounding non-zero histogram elements in             *
* histogram box.                                                           *
*                                                                          *
****************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::findBoundBoxForTags(
   hier::Box<DIM>& bound_box,
   const hier::IntVector<DIM>& min_box) const
{
   hier::Index<DIM> box_lo;
   hier::Index<DIM> box_hi;

   /*
    * Compute extent of bounding box in each coordinate direction.
    */
   for (int id = 0; id < DIM; id++) {
      boundTagHistogram(box_lo(id), box_hi(id), id, min_box(id));
   }

   bound_box = hier::Box<DIM>(box_lo, box_hi);
}

/*
****************************************************************************
*                                                                          *
* Bound the non-zero histogram elements within the interval                *
* (box_lo, box_hi) in the given coordinate direction. Note that            *
* computed bounding interval must be larger than given min_size.           *
*                                                                          *
* Note that it is assumed that box indices are cell indices.               *
*                                                                          *
****************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::boundTagHistogram(
   int& box_lo,
   int& box_hi,
   const int id,
   const int min_size) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (id >= 0) && (id < DIM) );
   for (int i = 0; i < id; i++) 
      assert(d_histogram[i].getSize() == d_box.numberCells(i));
#endif

   const int hist_lo = box_lo = d_box.lower(id);
   const int hist_hi = box_hi = d_box.upper(id);

   int width = box_hi-box_lo+1;  // Note: box indices are cell-centered

   if (width > min_size) {

      while ( (box_lo <= box_hi) &&
              (histogramElt(id, box_lo) == 0) ) box_lo++;
      while ( (box_hi >= box_lo) &&
              (histogramElt(id, box_hi) == 0) ) box_hi--;

      width = box_hi-box_lo+1;  // Note: box indices are cell-centered

      const int pad = min_size - width;
      if (pad > 0) {

         box_lo -= (pad+1)/2;
         if (box_lo < hist_lo) box_lo = hist_lo;

         box_hi = box_lo+min_size-1;
         if (box_hi > hist_hi) {
            box_hi = hist_hi;
            box_lo = tbox::Utilities::imax(hist_lo, box_hi-min_size+1);
         }

      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set the method used to perform historgram reduction.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::setCommunicationMode(int mode)
{
   d_communication_mode = mode;
}

/*
*************************************************************************
*                                                                       *
* Set the communicator to use during histogram reduction.               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::setCommunicator(tbox::MPI::comm comm)
{
   d_comm = comm;
}



/*
****************************************************************************
*                                                                          *
* Calculate histograms for tags matching specified value on patch level    *
* that lie within the bounds of the histogram box. Return the total        *
* number of tags within the box.                                           *
*                                                                          *
****************************************************************************
*/


template<int DIM> int HistogramBox<DIM>::computeTagHistogram(
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int tag_index,
   const int tag_val)
{
   if (d_communication_mode == ORIGINAL_MODE) {
     return(computeTagHistogramOrig(level, tag_index, tag_val));
   }

   t_compute->start();
   resetHistogram();
   for (int id = 0; id < DIM; id++) {
      computeLocalHistogram(id, level, tag_index, tag_val);
   }
   t_compute->stop();

   t_commwait->start();
   reduceTags(level);
   t_commwait->stop();

   int num_tags = 0;

   /*
    * the root processor counts the number of tags,
    * and broadcasts the number to the other processors.
    */
   const tbox::Array<int>& hist = histogram(0);
   const int hi = d_box.numberCells(0);
   for (int ic = 0; ic < hi; ic++) {
      num_tags += hist[ic];
   }

   if (d_communication_mode == BINARY_TREE_MODE) {
      t_commwait->start();
      tbox::Pointer< hier::BinaryTree<DIM> > tree = level->getBinaryTree();
      tree->partialBcast(d_box, num_tags);
      t_commwait->stop();
   }
   else if (d_communication_mode == COMMUNICATOR_MODE) {
      tbox::MPI::comm keep_me = tbox::MPI::getCommunicator();
      tbox::MPI::setCommunicator(d_comm);
      num_tags = tbox::MPI::bcast(num_tags, 0);
      tbox::MPI::setCommunicator(keep_me);
   }

   return(num_tags);
}



/*
*************************************************************************
*                                                                       *
* original version (prior to July 2002)                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> int HistogramBox<DIM>::computeTagHistogramOrig(
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int tag_index,
   const int tag_val)
{
   t_compute->start();
   resetHistogram();
   t_compute->stop();

   t_commwait->start();
   for (int id = 0; id < DIM; id++) {
      reduceTags(id, level, tag_index, tag_val);
   }
   t_commwait->stop();

   t_compute->start();

   int num_tags = 0;

   const tbox::Array<int>& hist = histogram(0);
   const int hi = d_box.numberCells(0);
   for (int ic = 0; ic < hi; ic++) {
      num_tags += hist[ic];
   }

   t_compute->stop();

   return(num_tags);
}




/*
****************************************************************************
*                                                                          *
* Calculate local histogram for given index space coordinate direction.    *
*                                                                          *
****************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::computeLocalHistogram(
   const int id,
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int tag_index,
   const int tag_val)
{
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      hier::Box<DIM> intersection = patch->getBox() * d_box;

      if ( !(intersection.empty()) ) {

         tbox::Pointer< pdat::CellData<DIM,int> >
            tag_data = patch->getPatchData(tag_index);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!tag_data.isNull());
#endif

         int ilo = intersection.lower(id);
         int ihi = intersection.upper(id);

         for (int ic_sb = ilo; ic_sb <= ihi; ic_sb++) {

            int tag_count = 0;

            hier::Box<DIM> src_box(intersection);
            src_box.lower(id) = ic_sb;
            src_box.upper(id) = ic_sb;
            for (pdat::CellIterator<DIM> ic(src_box); ic; ic++) {
               if ( (*(tag_data.getPointer()))(ic()) == tag_val ) {
                  tag_count++;
               }
            }

            histogramElt(id, ic_sb) += tag_count;
         }
      }
   }
}

/*
*****************************************************************************
*                                                                           *
* Calculate histogram for given index space coordinate direction.           *
* Note that this is a two step process. First, the local histogram          *
* values are computed for each patch that intersects the histogram box.     *
* Second, the local histogram arrays are reduced to determine the           *
* total histogram values across the entire patch level.                     *
* Note that the reduction is needed since the local histogram box           *
* will be contained within single patches but the global histogram          *
* box may span multiple patches that lie on different processors.           *
* The reduction will be performed on the tags lying within patch            *
* interiors only - ghost regions are not considered, to avoid multiple      *
* tag counting.                                                             *
*                                                                           *
*****************************************************************************
*/

template<int DIM> void HistogramBox<DIM>::reduceTags(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
   /*
    * Pack the entire histogram into a single buffer
    */
   int buf_size = 0;
   int id;
   for (id = 0; id < DIM; id++) {
      buf_size +=  d_box.numberCells(id);
   }

   tbox::Array<int> buffer;
   buffer.resizeArray(buf_size);
   int *buf_ptr = buffer.getPointer();

   int offset = 0;
   for (id = 0; id < DIM; id++) {
      int len = d_box.numberCells(id);
      memcpy(buf_ptr+offset, histogram(id).getPointer(), len*sizeof(int));
      offset += len;
   }

   /*
    * Perform sum reduction to form the global histogram
    */
   if (d_communication_mode == BINARY_TREE_MODE) {
      tbox::Pointer< hier::BinaryTree<DIM> > tree = level->getBinaryTree();
      tree->reduce(d_box, buf_ptr, buf_size); 
   }

   else if (d_communication_mode == COMMUNICATOR_MODE) {

      /*
       * perform an all-to-one reduction, but only amongst the processors
       * that are members of the communicator.
       */

      tbox::MPI::comm keep_me = tbox::MPI::getCommunicator();
      tbox::MPI::setCommunicator(d_comm);
      tbox::MPI::allToOneSumReduction(buf_ptr, buf_size);
      tbox::MPI::setCommunicator(keep_me);
   }

   else {  //d_communication_mode = ORIG_FAST_REDUCE
      tbox::MPI::sumReduction(buf_ptr, buf_size);
   }

   /*
    * Unpack the buffer
    */
   offset = 0;
   for (id = 0; id < DIM; id++) {


      int len = d_box.numberCells(id);
      memcpy(histogram(id).getPointer(), buf_ptr+offset, len*sizeof(int));
      offset += len;
   }
}

/* 
*************************************************************************
*                                                                       *
* This is the "original" version (prior to July 2002)                   *
*                                                                       *
*************************************************************************
*/
template<int DIM> void HistogramBox<DIM>::reduceTags(
   const int id,
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int tag_index,
   const int tag_val)
{
   t_compute->start();
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      hier::Box<DIM> intersection = patch->getBox() * d_box;

      if ( !(intersection.empty()) ) {
         tbox::Pointer< pdat::CellData<DIM,int> >
            tag_data = patch->getPatchData(tag_index);
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!tag_data.isNull());
#endif

         int ilo = intersection.lower(id);
         int ihi = intersection.upper(id);

         for (int ic_sb = ilo; ic_sb <= ihi; ic_sb++) {

            int tag_count = 0;

            hier::Box<DIM> src_box(intersection);
            src_box.lower(id) = ic_sb;
            src_box.upper(id) = ic_sb;
            for (pdat::CellIterator<DIM> ic(src_box); ic; ic++) {
               if ( (*(tag_data.getPointer()))(ic()) == tag_val ) {
                  tag_count++;
               }
            }

            histogramElt(id, ic_sb) += tag_count;
         }
      }
   }
   t_compute->stop();

   /*
    * After determining histogram for box in specified coordinate
    * direction on each local patch, add elements of histogram across
    * processors to obtain total for all patches on level.
    */

   t_commwait->start();
   tbox::MPI::sumReduction(histogram(id).getPointer(), 
                          d_box.numberCells(id));
   t_commwait->stop();

}




/*
***************************************************************************
*                                                                         *
* Release static timers.  To be called by shutdown registry to make sure  *
* memory for timers does not leak.                                        *
*                                                                         *
***************************************************************************
*/
template<int DIM>
void HistogramBox<DIM>::freeTimers()
{
   t_compute.setNull();
   t_commwait.setNull();
   return;
}

}
}
#endif
