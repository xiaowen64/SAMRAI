//
// File:        HistogramBox.h
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 296 $
// Modified:    $Date: 2005-04-14 16:39:17 -0700 (Thu, 14 Apr 2005) $
// Description: Histogram class for computing tagged cell signatures.
//

#ifndef included_mesh_HistogramBox
#define included_mesh_HistogramBox

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_MPI
#include "tbox/MPI.h"
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif


namespace SAMRAI {
    namespace mesh {

/*!
 * @brief Class HistogramBox<DIM> manages the histogram signature array for a
 * region of index space contained within a single box.  It is used during
 * the regridding process to gather information about the distribution of
 * tagged cells across a portion of a patch level 
 */

template<int DIM> class HistogramBox
{
public:
   /*!
    * Enumerated type for the communication mode.  For a description
    * of each of these modes, @see mesh::BergerRigoutsos.
    *
    * - \b ORIGINAL_MODE  {Original communication mode used in SAMRAI,
    *                      through v1.3.1 (global reductions).}
    * - \b ORIG_FAST_REDUCE_MODE {Original communication mode but with 
    *                             a 'faster' reduce operation.}
    * - \b BINARY_TREE_MODE {Reduces all-to-all communication costs
    *                        using explicit tbox::MPI send/recvs in a 
    *                        binary tree, rather than global reductions.}
    * - \b COMMUNICATOR_MODE {Similar to the binary tree option, but uses
    *                         tbox::MPI communicators in place of the explicit
    *                         send/recvs.}
    */
   enum{    ORIGINAL_MODE = 0,
            ORIG_FAST_REDUCE_MODE = 1, 
            BINARY_TREE_MODE = 2, 
            COMMUNICATOR_MODE = 3 };

   /*!
    * Default constructor for HistogramBox<DIM>. Note that the
    * HistogramBox<DIM> must be initialized with a box before
    * it can be used.
    */
   HistogramBox();

   /*!
    * Constructor for the HistogramBox<DIM>.  The HistogramBox<DIM>
    * will be initialized with the box and the histogram elements will
    * be allocated and initialized to zero.
    */
   HistogramBox(const hier::Box<DIM>& box);

   /**
    * Destructor for the HistogramBox.
    */
   ~HistogramBox<DIM>();

   /*!
    * Given histogram information in a box region, determine smallest
    * box containing non-zero histogram points.  Note that this routine
    * assumes that histogram values are already computed.
    */
   void findBoundBoxForTags(hier::Box<DIM>& bound_box,
                            const hier::IntVector<DIM>& min_box) const;

   /*!
    * Bound the non-zero histogram elements within the specified
    * interval (box_lo, box_hi) in the given coordinate direction.
    */
   void boundTagHistogram(int& box_lo,
                          int& box_hi,
                          const int id,
                          const int min_size) const;

   /*!
    * Set the method used to perform historgram reduction.
    * (This function should only by called from within
    * BergerRigoutsos<DIM>::findBoxesContainingTags.)
    */
   void setCommunicationMode(int mode);

   /*!
    * Set the communicator to use during histogram reduction.
    * This call only has an effect when used in conjunction with
    * setCommunicationMode(COMMUNICATOR_MODE);
    * (This function should only by called from within
    * BergerRigoutsos<DIM>::findBoxesContainingTags.)
    */
   void setCommunicator(tbox::MPI::comm comm);

   /*!
    * Compute histogram for integer tags (matching tag_val) within
    * the histogram box on the patch level.
    */
   int computeTagHistogram(const tbox::Pointer< hier::PatchLevel<DIM> > level,
                           const int tag_index,
                           const int tag_val);

//@{
   /*!
    * Compute global histogram by performing a reduction sum
    * on the local histograms.  Note that the
    * reduction is needed since the histogram box will not be contained
    * within any single patch, in general. Also, if the integer tag
    * patch data corresponding to the given descriptor index has ghost
    * cells, the values in these cells are ignored so that the reduction
    * sums data over patch interiors only.
    */
   void reduceTags(const tbox::Pointer< hier::PatchLevel<DIM> > level);

   void reduceTags(const int id, 
                   const tbox::Pointer< hier::PatchLevel<DIM> > level,
                   const int tag_index, 
                   const int tag_val);
//@}


   /*!
    * Return const reference to box represented by histogram.
    */
   const hier::Box<DIM>& getBox() const;
  
   /*!
    * Return non-const reference to histogram array for given direction.
    * Note that there is no bounds-checking on direction number.
    */
   tbox::Array<int>& histogram(const int id);

   /*!
    * Return const reference to histogram array for given direction.
    * Note that there is no bounds-checking on direction number.
    */
   const tbox::Array<int>& histogram(const int id) const;

   /*!
    * Return non-const histogram element corresponding to given
    * direction and cell index. Note cell index and histogram array
    * location are distinct, in general. Also, there is no bounds
    * checking done to assure the cell index lies within the box.
    */
   int& histogramElt(const int id, const int ic);

   /*!
    * Return const histogram element corresponding to given
    * direction and cell index. Note cell index and histogram array
    * location are distinct, in general. Also, there is no bounds
    * checking done to assure the cell index lies within the box.
    */
   const int& histogramElt(const int id, const int ic) const;

   /*!
    * Set the box for the HistogramBox, allocate storage for the histogram
    * signature arrays, and initialize the histogram elements to zero.
    */
   void setBox(const hier::Box<DIM>& box);

   /*!
    * Set all histogram elements to zero.
    */
   void resetHistogram();

private:

   int computeTagHistogramOrig(const tbox::Pointer< hier::PatchLevel<DIM> > level,
                               const int tag_index,
                               const int tag_val);

   void computeLocalHistogram(const int id,
                              const tbox::Pointer< hier::PatchLevel<DIM> > level,
                              const int tag_index,
                              const int tag_val);

   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();


   hier::Box<DIM> d_box;
   tbox::Array<int> d_histogram[DIM];

   tbox::MPI::comm d_comm;
   int d_communication_mode;

   static tbox::Pointer<tbox::Timer> t_compute;
   static tbox::Pointer<tbox::Timer> t_commwait;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "HistogramBox.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "HistogramBox.C"
#endif
