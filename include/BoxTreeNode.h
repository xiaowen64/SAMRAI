//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxTreeNode.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Basic building block used by BoxTree class.
//

#ifndef included_hier_BoxTreeNode
#define included_hier_BoxTreeNode

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
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


namespace SAMRAI {
   namespace hier {


/*!
 * @brief Building block used by BoxTree<DIM>.
 * 
 * This is a private class that is used by the BoxTree<DIM> class;
 * there is no reason that anyone should make direct use of this class.
 */

template<int DIM> class BoxTreeNode : public tbox::DescribedClass
{
public:

   struct Triple {
      Box<DIM> box;
      int       idx;
      int       owner;
   };

   /*!
    * @brief Constructs a BoxTreeNode that represents the physical
    * domain specified by the bounding box that includes all boxes
    * in \b box_triples.
    * 
    * Constructs a BoxTreeNode that represents the physical
    * domain specified by \b box.  The \b dim parameter is
    * the dimension along which the domain will be cut when
    * constructing child nodes.
    * 
    * @param box_triples input
    * @param have_mapping true is the calling BoxTree was passed a
    *                     processor mapping.
    * @param min_length if \b box_triples contains less than \b min_length
    *                   elements, than the list will never be partitioned
    *                   amongst child nodes.  Setting to a larger value
    *                   tends to decrease the total number of nodes in the tree
    *                   (and hence reduces memory requirements), but increase
    *                   the cost of findOverlappingBoxes.
    * @param dim the dimension along which the domain will be cut
    *            when constructing child nodes.
    * @param recurse_level is a counter that is used internally to
    *        keep track of the recurse level.  Do not change the
    *        default value.
    */
   BoxTreeNode(
      const tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples,
      bool have_mapping,
      int min_length = 10,
      int dim = DIM-1,
      int recurse_level = 0);

   /*!
    * The dtor does nothing even slightly interesting.
    */
   ~BoxTreeNode<DIM>();

   /*!
    * @brief Compute the box array indices of boxes that overlap the given \b box.
    *
    * If \b find_local_boxes = false, then the return array, \b indices,
    * contains the indices of the boxes (wrt the BoxArray that was passed
    * to the constructor) and that overlap with the specified \b box.
    * If \b find_local_boxes = true, then only the indices of those boxes
    * that overlap and whose corresponding patches are mapped to this
    * processor are returned.
    *
    * @param indices the indices of the overlapping boxes.
    * @param box the specified box whose overlaps are requested.
    * @param find_local_boxes switch to determine if all overlapping
    *        boxes are returned, or only those boxes that overlap
    *        and are local to this processor.
    * @param recurse_level is a counter that is used internally to
    *        keep track of the recurse level.  Do not change the
    *        default value.
    */
   void findOverlapIndices(tbox::List<int> &indices,
                           const Box<DIM> & box,
                           bool find_local_boxes,
                           int recurse_level = 0);

   /*!
    * @brief Create a list of all boxes that overlap the given \b box.
    *
    * @param overlap_boxes boxlist containing boxes that overlap with box.
    * @param box the specified box whose overlaps are requested.
    */
   void findOverlapBoxes(BoxList<DIM> &overlap_boxes,
                         const Box<DIM> & box);

   /*!
    * Counts the total number of nodes in all trees rooted at this
    * node, and add this to the totals storeed in BoxTree<DIM>::s_node_count
    * and BoxTree<DIM>::s_box_count.
    * This is a recursive call, i.e, all nodes in all subtrees
    * rooted at this node are also examined.
    * 
    * This call is primarily of interest for gathering information
    * on performance, during development and testing.
    */
   void count();

private:

   /*!
    * The copy ctor is not implemented.
    */
   BoxTreeNode(const BoxTreeNode<DIM>& box);

   /*!
    * The assignemnt operator is not implemented.
    */
   BoxTreeNode<DIM>& operator=(const BoxTreeNode<DIM>& box);

   /*!
    * The physical domain that this node represents.
    */
   Box<DIM> d_domain;

   /*!
    * Pointers to familial nodes.
    */
   tbox::Pointer< BoxTreeNode<DIM> > d_left_child;
   tbox::Pointer< BoxTreeNode<DIM> > d_right_child;

   /*!
    * The list of boxes, and their associated indices and processor
    * mappings, that are contained within the physical domain that 
    * this node represents.
    */
   tbox::Array<Triple> d_box_triples;

   /*!
    * Working space that is used in findOverlapIndices
    */
   tbox::Array<int> d_work;

   /*!
    * The dimension along which the input box triples will be partitioned.
    */
   int d_dim;

   /*!
    * A tree rooted at this node.  The tree is constructed
    * from d_box_triples; these boxes are partitioned along dimension
    * d_dim - 1.
    */
   tbox::Pointer< BoxTreeNode<DIM> > d_tree;

   /*!
    * This processor's rank.
    */
   int d_rank;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxTreeNode.C"
#endif

#endif

