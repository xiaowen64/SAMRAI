//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/clustering/BergerRigoutsos.h $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Berger/Rigoutsos tagged cell clustering algorithm.
//
 
#ifndef included_mesh_BergerRigoutsos
#define included_mesh_BergerRigoutsos

#include "SAMRAI_config.h"

#include "BoxGeneratorStrategy.h" 
#include "HistogramBox.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"

namespace SAMRAI {
    namespace mesh {

/*!
 * @brief BergerRigoutsos<DIM> provides operations that construct boxes
 * to cover a collection of tagged cells on a single AMR hierarchy patch
 * level.  This class is derived from the abstract base class
 * BoxGeneratorStrategy<DIM>.  Thus, it serves as a concrete implementation
 * of the box generator Strategy pattern interface.
 * 
 * The box generation algorithm is described in Berger and Rigoutsos,
 * tbox::IEEE Trans. on Sys, Man, and Cyber (21)5:1278-1286.
 *
 * NOTE: Algorithmic variations which may affect performance are available 
 * by calling the the static method BergerRigoutsos<DIM>::setClusteringOption(),
 * which sets the option for all instances of the class.
 *
 * @see mesh::HistogramBox
 * @see mesh::BoxGeneratorStrategy
 */

template<int DIM> class BergerRigoutsos 
: 
public BoxGeneratorStrategy<DIM>
{
public:
   /*!
    * Static function to set tag reduction option for clustering algorithm.
    *
    * @param  method   string identifying box intersection method.  Valid
    *                  choices are:  "ORIGINAL" (default case for single 
    *                  processor), "ORIG_FAST_REDUCE", "COMMUNICATOR", and
    *                  "BINARY_TREE" (default case for multiple processors).
    *                  The default empty string argument resets the method
    *                  to the default case.
    *
    * Each of these four options produces identical results but may result in
    * different performance characteristics depending on the machine architecture 
    * on which the code is run.  The defaults appear to be the most efficient 
    * choices for the tests we have performed.  The options are:
    *
    *   - ORIGINAL - original algorithm used in SAMRAI through v1.3.1;
    *                uses global reductions to collect tag information.
    *   - ORIG_FAST_REDUCE - original algorithm but with a faster reduce
    *                        operation.
    *   - COMMUNICATOR - uses communicators to store data in each recursion
    *                    of the algorithm.
    *   - BINARY_TREE - uses a binary tree communication scheme with
    *                   hand-coded sends and receives.
    * 
    * If an invalid non-empty string is passed, an unrecoverable error will result.
    */
   static void setClusteringOption(const std::string& method = std::string());

   /*!
    * Virtual destructor.
    */
   virtual ~BergerRigoutsos<DIM>();


   /*!
    * The constructor queries the tbox::InputManager to determine
    * which of the various findBoxesContainingTags algorithms
    * should be used.
    */
   BergerRigoutsos();

   /*!
    * Create a list of boxes that covers all integer tags on the patch level
    * that match the specified tag value.  Each box will be at least as
    * large as the given minimum size and the tolerances will be met.  The
    * 
    * The efficiency tolerance is a threshold value for the percentage of
    * tagged cells in each box.  If this percentage is below the tolerance,
    * the box will continue to be split into smaller boxes.
    * 
    * The combine tolerance is a threshold value for the sum of the volumes
    * of two boxes into which a box may be potentially split.  If ratio of
    * that sum and the volume of the original box, the box will not be split.
    */
   void findBoxesContainingTags(hier::BoxList<DIM>& boxes,
                                const tbox::Pointer< hier::PatchLevel<DIM> > level,
                                const int index,
                                const int tag_val,
                                const hier::Box<DIM>& bound_box,
                                const hier::IntVector<DIM>& min_box,
                                const double efficiency_tol,
                                const double combine_tol) const; 


private:
   /*
    * This algorithm invokes the "original" implementation 
    * used in SAMRAI through v1.3.1.  It uses global reductions
    * to collect tag information.  It is possible to use a faster
    * implementation of global reductions by choosing the "fast_reduce"
    * option.  This method will be invoked if the user calls 
    * setClusteringOption("ORIGINAL") or setClusteringOption("ORIG_FAST_REDUCE").
    * This is the default algorithm choice in the single processor case.
    */
   void findBoxesContainingTagsOriginal(
      hier::BoxList<DIM>& boxes,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const int tag_val,
      const hier::Box<DIM>& bound_box,
      const hier::IntVector<DIM>& min_box,
      const double efficiency_tol,
      const double combine_tol) const; 

   /*
    * This algorithm reduces all-to-all communication costs by forming
    * communicators on recursive calls.  Only processors that have local
    * patches that intersect with the current bound_box are members
    * of the new communicator and participate in the recursion.
    * This method will be invoked if the user calls
    * setClusteringOption("COMMUNICATOR").
    */
   void findBoxesContainingTagsCommunicators(
      hier::BoxList<DIM>& boxes,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const int tag_val,
      const hier::Box<DIM>& bound_box,
      const hier::IntVector<DIM>& min_box,
      const double efficiency_tol,
      const double combine_tol,
      int recurse_level,
      tbox::SAMRAI_MPI::comm comm) const;

   /*
    * This method reduces all-to-all communication costs in a fashion
    * similar to the previous method; here, however, all communications
    * are explicitly managed using point-to-point communications and
    * a binary tree, instead of forming new communicators. This method 
    * will be invoked if the user calls setClusteringOption("BINARY_TREE").
    * This is the default in the case of more than one processor.
    */
   void findBoxesContainingTagsBinaryTree(
      hier::BoxList<DIM>& boxes,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const int tag_val,
      const hier::Box<DIM>& bound_box,
      const hier::IntVector<DIM>& min_box,
      const double efficiency_tol,
      const double combine_tol,
      int recurse_level) const;

   /*
    * Attempt to split the given box which bounds a collection of tagged
    * cells into two smaller boxes, eliminating non-tagged cells from the
    * union.  If an appropriate splitting is found, the two smaller boxes
    * boxes are returned and the return value of the function is true.
    * Otherwise, false is returned and the argument boxes remain empty.
    */
   bool splitTagBoundBox(hier::Box<DIM>& box_lft,
                         hier::Box<DIM>& box_rgt,
                         const hier::Box<DIM>& bound_box,
                         const HistogramBox<DIM>& hist_box,
                         const hier::IntVector<DIM>& min_box) const;

   /*
    * Attempt to find a zero value in the tag histogram near the middle of
    * the index interval (lo, hi) in the given coordinate direction.  If
    * such a point is found, true is returned; otherwise, false is returned.
    */
   bool findZeroCutPoint(int& cut_pt,
                         const int id,
                         const int lo,
                         const int hi,
                         const HistogramBox<DIM>& hist_box,
                         const int min_size) const;

   /*
    * Attempt to find a point in the given coordinate direction near an
    * inflection point in the tag histogram for that direction.  The search
    * is restricted to the sub-interval (lo, hi).  If no such inflection
    * point is found, the mid-point of the interval is returned as the
    * cut point.
    */
   void cutAtLaplacian(int& cut_pt,
                       const int id,
                       const int lo,
                       const int hi,
                       const HistogramBox<DIM>& hist_box,
                       const int min_size) const;

   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();

   /*
    * Choice of algorithm for tag reductions.  
    */
   static std::string s_tag_reduce_method;
   static bool s_tag_reduce_method_initialized;

   // The following are not implemented:
   BergerRigoutsos(const BergerRigoutsos<DIM>&);
   void operator=(const BergerRigoutsos<DIM>&);

   static tbox::Pointer<tbox::Timer> t_cluster;
   static tbox::Pointer<tbox::Timer> t_globalize_boxes;
   static tbox::Pointer<tbox::Timer> t_compute;
   static tbox::Pointer<tbox::Timer> t_commwait;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BergerRigoutsos.C"
#endif
