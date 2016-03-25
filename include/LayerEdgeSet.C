/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/dlbg/LayerEdgeSet.C $
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1846 $
 * Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
 * Description: Set of edges in distributed box graph.
 */

#ifndef included_hier_LayerEdgeSet_C
#define included_hier_LayerEdgeSet_C

#include "LayerEdgeSet.h"
#include "BoxTree.h"

#include "tbox/SAMRAI_MPI.h"

#include IOSTREAM_HEADER_FILE
#include IOMANIP_HEADER_FILE


#ifdef DEBUG_NO_INLINE
#include "LayerEdgeSet.I"
#endif

namespace SAMRAI {
namespace hier {


#define SIZE_EXCHANGE_TAG (100000)
#define EDGE_EXCHANGE_TAG (100001)

#define DEBUG_SEND_EDGE_EXTRA_DATA


template<int DIM>
LayerEdgeSet<DIM>::LayerEdgeSet()
   : d_gcw(0),
     d_head_ratio(0) /* nonsense value */,
     d_base(NULL),
     d_cnect(new Connectivity[1]),
     d_parallel_state(DISTRIBUTED),
     d_partner(NULL),
     d_partner_is_managed(false),
     d_rank(tbox::SAMRAI_MPI::getRank()),
     d_nprocs(tbox::SAMRAI_MPI::getNodes())
{
   TBOX_ASSERT( d_cnect != NULL );
   return;
}


template<int DIM>
LayerEdgeSet<DIM>::~LayerEdgeSet()
{
   deallocateDataSelf();
   d_base = NULL;
   delete [] d_cnect;
   d_cnect = NULL;
   if ( d_partner_is_managed ) {
      delete d_partner;
   }
   d_partner = NULL;
   return;
}


template<int DIM>
const typename LayerEdgeSet<DIM>::Connectivity
&LayerEdgeSet<DIM>::getConnectivity( const int rank ) const
{
   if ( d_parallel_state == DISTRIBUTED && rank > -1 && rank != d_rank ) {
      TBOX_ERROR("Non-local connectivity unavailable in distributed state.");
   }
   return d_parallel_state == DISTRIBUTED ? *d_cnect : d_cnect[rank>-1?rank:d_rank];
}


template<int DIM>
typename LayerEdgeSet<DIM>::Connectivity
&LayerEdgeSet<DIM>::getConnectivity( const int rank )
{
   if ( d_parallel_state == DISTRIBUTED && rank > -1 && rank != d_rank ) {
      TBOX_ERROR("Non-local connectivity unavailable in distributed state.");
   }
   return d_parallel_state == DISTRIBUTED ? *d_cnect : d_cnect[rank>-1?rank:d_rank];
}


template<int DIM>
void LayerEdgeSet<DIM>::setMaxGhostCellWidth(
   const hier::IntVector<DIM> &gcw )
{
   if ( d_gcw != gcw ) {
      deallocateData();
      d_gcw = gcw;
   }
   return;
}


template<int DIM>
LayerEdgeSet<DIM> &LayerEdgeSet<DIM>::getPartner()
{
  if ( d_partner == NULL ) {
    TBOX_ERROR("No partner");
  }
  return *d_partner;
}


template<int DIM>
const LayerEdgeSet<DIM> &LayerEdgeSet<DIM>::getPartner() const
{
  if ( d_partner == NULL ) {
    TBOX_ERROR("No partner");
  }
  return *d_partner;
}


template<int DIM>
void LayerEdgeSet<DIM>::attachPartner( LayerEdgeSet &partner )
{
  if ( d_partner != NULL || partner.d_partner != NULL ) {
    TBOX_ERROR("Cannot attach multiple partners");
  }
  if ( partner.d_parallel_state != d_parallel_state ) {
     TBOX_ERROR("Partners do not have same parallel states.");
  }
  if ( partner.getBaseRefinementRatio() != getHeadRefinementRatio() ||
       partner.getHeadRefinementRatio() != getBaseRefinementRatio() ) {
     TBOX_ERROR("Partners do not agree on head-tail refinement ratios.");
  }
  if ( (partner.d_gcw*getBaseRefinementRatio()) !=
       (d_gcw*partner.getBaseRefinementRatio()) ) {
     TBOX_ERROR("Partners do not have compatible overlap criteria.");
  }
  d_partner = &partner;
  d_partner_is_managed = false;
  partner.d_partner = this;
  partner.d_partner_is_managed = false;
  return;
}


template<int DIM>
LayerEdgeSet<DIM> &LayerEdgeSet<DIM>::createPartner(
   const LayerNodeSet<DIM> &partner_base,
   const Connectivity *partner_connectivity )
{
  if ( d_partner != NULL ) {
    TBOX_ERROR("Cannot create partner while another is attached");
  }
  if ( partner_base.getRefinementRatio() != d_head_ratio ) {
     TBOX_ERROR("LayerEdgeSet::createPartner(): partner base refinement ratio\n"
                <<"must be the same as the head refinement ratio.\n");
  }

  d_partner = new LayerEdgeSet<DIM>;
  TBOX_ASSERT( d_partner != NULL );

  hier::IntVector<DIM> partner_gcw =
     ( d_gcw * getBaseRefinementRatio() ) / getHeadRefinementRatio();
  TBOX_ASSERT( ( partner_gcw * getHeadRefinementRatio() ) ==
               ( d_gcw * getBaseRefinementRatio() ) );

  d_partner->initialize( DISTRIBUTED,
                         partner_base,
                         d_head_ratio,
                         partner_gcw,
                         partner_connectivity );

  d_partner_is_managed = true;
  d_partner->d_partner = this;
  d_partner->d_partner_is_managed = false;

  return *d_partner;
}


template<int DIM>
LayerEdgeSet<DIM> &LayerEdgeSet<DIM>::operator=( const LayerEdgeSet &r )
{
   d_gcw = r.d_gcw;
   d_base = r.d_base;
   d_parallel_state = d_parallel_state;

   if ( d_parallel_state == DISTRIBUTED ) {
      d_cnect = new Connectivity[1];
      *d_cnect = *r.d_cnect;
   }
   else {
      d_cnect = new Connectivity[d_nprocs];
      int n;
      for ( n=0; n<d_nprocs; ++n ) {
         d_cnect[n] = r.d_cnect[n];
      }
   }

   if ( r.d_partner == &r ) {
      attachPartner(*this);
   }
   else if ( r.d_partner != NULL ) {
      createPartner( *r.d_partner->d_base,
                     &r.d_partner->getConnectivity() );
      /*
       * To set data for a created partner,
       * do the assignment without touching its partner data,
       * which was just set up by the call to createPartner().
       */
      // d_partner->d_gcw = r.d_partner->d_gcw;
      // d_partner->d_base = r.d_partner->d_base;
      // d_partner->d_parallel_state = d_partner->d_parallel_state;
      // if ( d_partner->d_parallel_state == DISTRIBUTED ) {
         // d_partner->d_cnect = new Connectivity[1];
         // *d_partner->d_cnect = *r.d_partner->d_cnect;
      // }
      // else {
         // d_partner->d_cnect = new Connectivity[d_nprocs];
         // int n;
         // for ( n=0; n<d_partner->d_nprocs; ++n )  {
            // d_partner->d_cnect[n] = r.d_partner->d_cnect[n];
         // }
      // }
   }
   else {
      d_partner = NULL;
   }

   return *this;
}


template<int DIM>
void LayerEdgeSet<DIM>::deallocateData()
{
   deallocateDataSelf();
   if ( d_partner ) {
      d_partner->deallocateDataSelf();
   }
   return;
}


template<int DIM>
void LayerEdgeSet<DIM>::deallocateDataSelf()
{
   int n;
   switch (d_parallel_state) {
   case DISTRIBUTED:
      d_cnect->clear();
      break;
   case GLOBALIZED:
      for ( n=0; n<d_nprocs; ++n ) {
         d_cnect[n].clear();
      }
      break;
   default:
      TBOX_ERROR("Corrupted parallel state.\n");
   }
   return;
}


template<int DIM>
void LayerEdgeSet<DIM>::allocateConnectivity( const NodeContainer &nodes,
                                               Connectivity &cnect ) const
{
   typename NodeContainer::const_iterator i_node;
   NodeContainer empty_container;
   for ( i_node=nodes.begin(); i_node!=nodes.end(); ++i_node ) {
      cnect.insert( cnect.end(),
#ifdef USE_SORTED_SET
                    typename Connectivity::value_type((*i_node).getLocalIndex(),
                                             empty_container)
#else
                    typename Connectivity::value_type((*i_node).getLocalIndex(),
                                             empty_container)
#endif
                    );
   }
   return;
}





template<int DIM>
void LayerEdgeSet<DIM>::findEdges()
{
   if ( ! d_partner ) {
      TBOX_ERROR("No partner attached in findEdges()."
                 << "A partner is required before findEdges() can be called.");
   }
   findEdges( *this, *d_partner );
   return;
}

template<int DIM>
void LayerEdgeSet<DIM>::findEdges( LayerEdgeSet &a, LayerEdgeSet &b )
{
   if ( (a.d_gcw*b.getBaseRefinementRatio()) !=
        (b.d_gcw*a.getBaseRefinementRatio()) ) {
      TBOX_ERROR("Cannot guarantee consistent edge data without\n"
                 <<"consistent max ghost width.\n");
   }

   /*
    * We will use findEdges_...() to find the edges.
    * Layer nodes for both a and b must be GLOBALIZED.
    * If they are not, create temporary GLOBALIZED versions.
    */

   LayerNodeSet<DIM> *temp_globalized_base = NULL;
   LayerNodeSet<DIM> *temp_globalized_head = NULL;

   const LayerNodeSet<DIM> *orig_node_set_a = a.d_base;
   if ( a.d_base->getParallelState() == LayerNodeSet<DIM>::DISTRIBUTED ) {
      temp_globalized_base = copyAndGlobalize(*a.d_base);
      temp_globalized_base->setParallelState(LayerNodeSet<DIM>::GLOBALIZED);
      a.d_base = temp_globalized_base;
   }
   const LayerNodeSet<DIM> *orig_node_set_b = b.d_base;
   if ( &a != &b &&
        b.d_base->getParallelState() == LayerNodeSet<DIM>::DISTRIBUTED ) {
      temp_globalized_head = copyAndGlobalize(*b.d_base);
      temp_globalized_head->setParallelState(LayerNodeSet<DIM>::GLOBALIZED);
      b.d_base = temp_globalized_head;
   }

   /*
    * findEdges is not symmetric between a and b,
    * so call it for both..
    * One exception is if b is actually *this, in which case,
    * symmetry is automatic.
    */
   if(1) a.findEdges_rbbt( *b.d_base, a.d_base == b.d_base );
   else a.findEdges_quadratic( b );
#ifdef DEBUG_CHECK_ASSERTIONS
   a.checkNodeNabrCorrespondance();
#endif
   if ( &b != &a ) {
      if (1) b.findEdges_rbbt( *a.d_base, a.d_base == b.d_base );
      else b.findEdges_quadratic( a );
#ifdef DEBUG_CHECK_ASSERTIONS
      b.checkNodeNabrCorrespondance();
#endif
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   a.checkNodeNabrCorrespondance();
   b.checkNodeNabrCorrespondance();
#endif

   /*
    * Restore original layer node sets and deallocate temporary ones.
    */
   a.d_base = orig_node_set_a;
   if ( &a != &b ) {
      b.d_base = orig_node_set_b;
   }
   if ( temp_globalized_base != NULL ) {
      delete temp_globalized_base;
   }
   if ( temp_globalized_head != NULL ) {
      delete temp_globalized_head;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   a.checkNodeNabrCorrespondance();
   b.checkNodeNabrCorrespondance();
#endif

   return;
}




template<int DIM>
void LayerEdgeSet<DIM>::findEdges_rbbt( const LayerNodeSet<DIM> &head,
                                        const bool ignore_self_overlap )
{
   /*
    * Finding edges for this object, using
    * an externally provided head LayerNodeSet.
    *
    * Global nodes provided by other are placed in a hier::BoxTree<DIM>
    * so they can be quickly searched to see which intersects the
    * boxes in this object.
    */
   if ( head.getParallelState() != LayerNodeSet<DIM>::GLOBALIZED ) {
      TBOX_ERROR("LayerEdgeSet::findEdges_rbbt() requires given head\n"
                 <<"to be GLOBALIZED.\n");
   }


   /*
    * Explicitly generate the one-to-one mapping between the nodes
    * and their neighbor lists.  As we find edges in this method,
    * the mappings are implicitly generated only for nodes that have
    * neighbors.
    */
   allocateConnectivity( getNodeContainer(), getConnectivity() );

   /*
    * The nomenclature "sthis" refers to the *this layer
    * and "other" refer to the layer in the argument.
    */

   /*
    * Determine relationship between base and head index spaces.
    */
   const bool other_is_finer =
      head.getRefinementRatio() > d_base->getRefinementRatio();
   const hier::IntVector<DIM> other_sthis_ratio = other_is_finer ?
      head.getRefinementRatio()/d_base->getRefinementRatio() :
      hier::IntVector<DIM>(0);

   const bool sthis_is_finer =
      d_base->getRefinementRatio() > head.getRefinementRatio();
   const hier::IntVector<DIM> sthis_other_ratio = sthis_is_finer ?
      d_base->getRefinementRatio()/head.getRefinementRatio() :
      hier::IntVector<DIM>(0);

   typename NodeContainer::const_iterator node_iterator;

   /*
    * Create a box array version of the head layer
    * so we can use hier::BoxTree<DIM>.
    * This implicitly creates a global index for each of the global node.
    * Create a map from the global index back to the nodes in the node
    * container so we can map the hier::BoxTree result back to layer nodes.
    */
   int box_array_size = 0;
   int n;
   for ( n=0; n<d_nprocs; ++n ) {
      const NodeContainer &other_node_container = head.getNodeContainer(n);
      box_array_size += other_node_container.size();
   }
   hier::BoxArray<DIM> box_array(box_array_size);
   Node const **index_to_node = new Node const *[box_array_size];
   int i = 0;
   for ( n=0; n<d_nprocs; ++n ) {
      const NodeContainer &other_node_container = head.getNodeContainer(n);
      for ( node_iterator=other_node_container.begin();
            node_iterator!=other_node_container.end();
            ++node_iterator ) {
         const Node &node = *node_iterator;
         index_to_node[i] = &node;
         box_array(i) = node.getBox();
         if ( sthis_is_finer ) {
            box_array(i).refine(sthis_other_ratio);
         }
         else {
            hier::IntVector<DIM> equiv_head_gcw =
               d_gcw*head.getRefinementRatio()/getBaseRefinementRatio();
            box_array(i).grow(equiv_head_gcw);
         }
         ++i;
      }
   }
   hier::BoxTree<DIM> box_tree(box_array);


   /*
     A neighbor of a node would be discarded if
     - the two are equal by comparison,
     - they are from layers with the same refinement ratio, and
     - ignore_self_overlap is true.
   */
   const bool discard_self_overlap =
      ignore_self_overlap &&
      ( d_base->getRefinementRatio() == head.getRefinementRatio() );


   /*
    * Use hier::BoxTree<DIM> to find intersecting boxes in the other layer.
    * Use index_to_node to map the indices back to nodes and save data.
    */
   const NodeContainer &node_container = getNodeContainer();
   Connectivity &node_connectivity = getConnectivity();
   tbox::Array<int> overlap_indices;
   for ( node_iterator=node_container.begin();
         node_iterator!=node_container.end();
         ++node_iterator ) {

      const Node &node = *node_iterator;
      NabrContainer &nabrs_for_box = node_connectivity[node.getLocalIndex()];

      hier::Box<DIM> box = node.getBox();
      if ( other_is_finer ) box.refine(other_sthis_ratio);
      else if ( sthis_is_finer ) box.grow(d_gcw);
      box_tree.findOverlapIndices( overlap_indices, box );

      const int num_overlaps = overlap_indices.size();
      nabrs_for_box.clear();
#ifdef USE_SORTED_SET
      nabrs_for_box.reserve(num_overlaps);
#endif
      for ( i=0; i<num_overlaps; ++i ) {

         const Node &nabr = *(index_to_node[overlap_indices[i]]);

         /*
           By convention, we do not make a node its own neighbor.
           node and nabr are the same node only if they are the
           same by comparison AND they are part of the same set.
         */
         if ( discard_self_overlap && node == nabr ) continue;

         /*
          * The box_array was sorted in the correct node ordering convention,
          * so the box_tree should return the results in the correct ordering
          * also.  Therefore, we can append the nodes sequentially (without
          * searching for their positions in the neighbor container.
          */
#ifdef USE_SORTED_SET
         nabrs_for_box.append(nabr);
#else
         nabrs_for_box.insert( nabrs_for_box.end(), nabr );
#endif
      }
   }

   delete [] index_to_node;

   return;
}




#if 1
template<int DIM>
void LayerEdgeSet<DIM>::findEdges_quadratic( const LayerEdgeSet &other )
{
   /*
    * Finding edges for this object, using
    * an another layeredge set whose layer node set is in hybrid mode.
    * Locally examining the entire layer,
    * including local and non-local nodes in the layer.
    */
   if ( d_base->getParallelState() != LayerNodeSet<DIM>::GLOBALIZED ) {
      TBOX_ERROR("LayerEdgeSet::findEdges() can only be used\n"
                 <<"when parallel state is Duplicated.\n"
                 <<"Use LayerNodeSet::setParallelState(GLOBALIZED) to do this.\n");
   }


   /*
    * Explicitly generate the one-to-one mapping between the nodes
    * and their neighbor lists.  As we find edges in this method,
    * the mappings are implicitly generated only for nodes that have
    * neighbors.
    */
   allocateConnectivity( getNodeContainer(), getConnectivity() );

   /*
    * The nomenclature "sthis" refers to the *this layer
    * and "other" refer to the layer in the argument.
    */

   const NodeContainer &_sthis_nodes = getNodeContainer();
   Connectivity &sthis_cnect = this->getConnectivity();
   hier::IntVector<DIM> sthis_gcw = this->d_gcw;

   const bool other_is_finer =
      other.getBaseRefinementRatio() > this->getBaseRefinementRatio();
   const hier::IntVector<DIM> other_sthis_ratio = other_is_finer ?
      other.getBaseRefinementRatio()/this->getBaseRefinementRatio() : hier::IntVector<DIM>(0);

   const bool sthis_is_finer =
      this->getBaseRefinementRatio() > other.getBaseRefinementRatio();
   const hier::IntVector<DIM> sthis_other_ratio = sthis_is_finer ?
      this->getBaseRefinementRatio()/other.getBaseRefinementRatio() : hier::IntVector<DIM>(0);


   /*
    * Create a copy of the list of neighbors to the sthis-loop layer,
    * so we can modify (either refine or grow) its boxes for the purpose
    * of intersection check.
    * For intersection checks, the boxes must be compared at
    * the same index space and the node on the finer layer
    * must be grown before checking.
    *
    * Note: it probably suffices to dupplicate just the boxes for
    * modification, not the whole node.  For now, copy the whole node.
    */
   NabrContainer sthis_nodes(_sthis_nodes);
   hier::BoxArray<DIM> sthis_boxes(_sthis_nodes.size());
   int seq_ndx;
   if ( other_is_finer ) {
      typename NodeContainer::iterator j;
      for ( seq_ndx=0, j=sthis_nodes.begin();
            j!=sthis_nodes.end(); ++j, ++seq_ndx ) {
         sthis_boxes(seq_ndx) = (*j).getBox();
         sthis_boxes(seq_ndx).refine(other_sthis_ratio);
         // (*j).getBox().refine(other_sthis_ratio);
      }
   }
   else {
      typename NabrContainer::iterator j;
      for ( seq_ndx=0, j=sthis_nodes.begin();
            j!=sthis_nodes.end(); ++j, ++seq_ndx ) {
         sthis_boxes(seq_ndx) = (*j).getBox();
         sthis_boxes(seq_ndx).grow(d_gcw);
         // (*j).getBox().grow(d_gcw);
      }
   }

   /*
    * Loop through all nodes in the other layer, searching for
    * intersection with the sthis layer.
    */

   int n;
   typename NodeContainer::const_iterator i_other_node;
   for ( n=0; n<d_nprocs; ++n ) {
      for ( i_other_node=other.getNodeContainer(n).begin();
            i_other_node!=other.getNodeContainer(n).end();
            ++i_other_node ) {

         const Node &other_node = *i_other_node;

         hier::Box<DIM> other_box = other_node.getBox();
         if ( other_is_finer ) other_box.grow(d_gcw);
         else if ( sthis_is_finer ) other_box.refine(sthis_other_ratio);

         typename NodeContainer::const_iterator i_sthis_node;
         for ( seq_ndx=0, i_sthis_node=sthis_nodes.begin();
               i_sthis_node!=sthis_nodes.end(); ++seq_ndx, ++i_sthis_node ) {

            const Node &sthis_node = *i_sthis_node;
            // hier::Box<DIM> sthis_box = sthis_node.getBox();
            const hier::Box<DIM> &sthis_box = sthis_boxes(seq_ndx);

            if ( other_box.intersects(sthis_box) ) {
               /*
                * sthis_node_nabrs is the neighbor container for sthis_node,
                * if sthis_node is a local node.  We do not set neighbor info
                * for non-local nodes.
                */
               NabrContainer &sthis_node_nabrs =
                  sthis_cnect[sthis_node.getLocalIndex()];

#ifdef USE_SORTED_SET
               sthis_node_nabrs.append( other_node );
#else
               sthis_node_nabrs.insert( sthis_node_nabrs.end(), other_node );
#endif
            }

         }

      }
   }

   return;
}
#endif




template<int DIM>
void LayerEdgeSet<DIM>::findEdges_serial( LayerEdgeSet &other )
{
   /*
    * Finding edges with no extra info, requires locally examining
    * the entire layer.
    */
   TBOX_ERROR("LayerEdgeSet::findEdges_serial() not yet rewriten\n"
              "since regrouping!  Not to mention totally untested!\n");
   if ( d_parallel_state != GLOBALIZED ) {
      TBOX_ERROR("LayerEdgeSet::findEdges_serial() can only be used\n"
                 <<"when the node parallel state is GLOBALIZED.\n"
                 <<"Use LayerNodeSet::setParallelState(GLOBALIZED) to\n"
                 <<"do this.\n");
   }


   /*
    * We treat *this and other symmetrically--it matters not if
    * they are switched.  But for efficiency, we would like to put
    * the outer search loop over the layer with the larger number
    * of nodes.
    *
    * We make a (deep) copy of the node container for the smaller
    * layer and make modifications to those nodes.  This eliminates
    * making the changes repeadedly with each new iteration of the
    * outer loop.
    *
    * The nomenclature "outer" and "inner" refer to the layer of
    * the outer or inner search loop.
    */

   LayerEdgeSet &outer_layer =
      getNodeContainer().size() > other.getNodeContainer().size() ? *this : other;
   LayerEdgeSet &inner_layer =
      getNodeContainer().size() <= other.getNodeContainer().size() ? *this : other;

   const NodeContainer &outer_nodes = outer_layer.getNodeContainer();
   Connectivity *outer_cnect = outer_layer.d_cnect;
   hier::IntVector<DIM> outer_gcw = outer_layer.d_gcw;

   const NodeContainer &_inner_nodes = inner_layer.getNodeContainer();
   Connectivity *inner_cnect = inner_layer.d_cnect;
   hier::IntVector<DIM> inner_gcw = inner_layer.d_gcw;

   const bool outer_is_finer =
      outer_layer.getBaseRefinementRatio() > inner_layer.getBaseRefinementRatio();
//   const hier::IntVector<DIM> outer_inner_ratio = outer_is_finer ?
//     outer_layer.getBaseRefinementRatio()/inner_layer.getBaseRefinementRatio() : hier::IntVector<DIM>(0);

   const bool inner_is_finer =
      inner_layer.getBaseRefinementRatio() > outer_layer.getBaseRefinementRatio();
   const hier::IntVector<DIM> inner_outer_ratio = inner_is_finer ?
     inner_layer.getBaseRefinementRatio()/outer_layer.getBaseRefinementRatio() : hier::IntVector<DIM>(0);

   /*
    * Create a copy of the list of neighbors to the inner-loop layer,
    * so we can modify it for the purpose of comparison.
    * We either refine these nodes or grow them.  For
    * intersection checks, the boxes must be compared at
    * the same index space and the node on the finer layer
    * must be grown before checking.
    *
    * Note: it probably suffices to dupplicate just the boxes for
    * modification, not the whole node.
    */
   NabrContainer inner_nodes(_inner_nodes);
   if ( outer_is_finer ) {
      typename NabrContainer::iterator j;
      for ( j=inner_nodes.begin();
            j!=inner_nodes.end(); ++j ) {
         // (*j).getBox().refine(outer_inner_ratio);
      }
   }
   else {
      typename NabrContainer::iterator j;
      for ( j=inner_nodes.begin();
            j!=inner_nodes.end(); ++j ) {
         // (*j).getBox().grow(d_gcw);
      }
   }

   typename NodeContainer::const_iterator i_outer_node;
   for ( i_outer_node=outer_nodes.begin();
         i_outer_node!=outer_nodes.end(); ++i_outer_node ) {

      const Node &outer_node = *i_outer_node;

      hier::Box<DIM> outer_box = outer_node.getBox();
      if ( outer_is_finer ) outer_box.grow(d_gcw);
      else if ( inner_is_finer ) outer_box.refine(inner_outer_ratio);

      const int outer_owner = outer_node.getOwnerRank();

      /*
       * outer_node_nabrs is the neighbor container for outer_node,
       * if outer_node is a local node.  We do not set neighbor info
       * for non-local nodes.
       */
      NabrContainer &outer_node_nabrs =
         outer_cnect[outer_node.getOwnerRank()][outer_node.getLocalIndex()];

      typename NodeContainer::const_iterator i_inner_node, i__inner_node;
      for ( i_inner_node=inner_nodes.begin(),i__inner_node=_inner_nodes.begin();
            i_inner_node!=inner_nodes.end(); ++i_inner_node,++i__inner_node ) {

         const Node &inner_node = *i_inner_node;
         hier::Box<DIM> inner_box = inner_node.getBox();

         const int inner_owner = inner_node.getOwnerRank();

         if ( (outer_owner == d_rank || inner_owner == d_rank)
              && outer_box.intersects(inner_box) ) {

            /*
             * inner_node_nabrs is the neighbor container for inner_node,
             * if inner_node is a local node.  We do not set neighbor info
             * for non-local nodes.
             */
            NabrContainer &inner_node_nabrs =
               inner_cnect[inner_node.getOwnerRank()][inner_node.getLocalIndex()];

            outer_node_nabrs.insert( *i__inner_node );
            inner_node_nabrs.insert( outer_node );
         }

      }

   }

   return;
}




/*
***********************************************************************
*                                                                     *
*                        +--------> (head layer)                      *
*                        |                                            *
*                        |                ^                           *
*         leg to head -> |                |                           *
*                        |                |                           *
*                                         |                           *
*                  (middle layer)         | <- bridge                 *
*                                         |                           *
*                        |                |                           *
*         leg to tail -> |                |                           *
*                        |                |                           *
*                        |                                            *
*                        +--------> (tail layer)                      *
*                                                                     *
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::bridge( const LayerEdgeSet &edge_to_head,
                                const LayerEdgeSet &edge_to_base )
{

   if ( edge_to_head.d_partner == NULL || edge_to_base.d_partner == NULL ) {
      TBOX_ERROR("Legs to head and tail must have partners in\n"
                 <<"LayerEdgeSet::bridge.\n");
   }


   /*
    * As the edges between the tail and head layers are discovered,
    * place the info directly into message buffers to be sent to
    * processes that need the info.
    */

   // map< int, EdgeSharingTransaction > mesg_to_owners;
   // bridge_discoverEdges( first_leg, secnd_leg, mesg_to_owners );
   // bridge_shareEdges( mesg_to_owners );

   std::map< int, CommunicationStruct > mesg_to_owners;

   if ( d_parallel_state != GLOBALIZED ) {

      /*
       * Create (empty) messages to owners of nodes in the head and tail
       * layers (except for local process).
       * This forces a message exchange with those processes in
       * bridge_shareEdges call.  This exchange is necessary, because
       * even if no edges are locally discovered for some nodes in the head
       * and tail layers, those owners may have discovered some edges to
       * share with the local process.
       *
       * This is not done in serial mode because no communication is needed.
       */
      std::set<int> owner_set;
      getNabrOwnerSet( edge_to_head.getConnectivity(), owner_set );
      getNabrOwnerSet( edge_to_head.d_partner->getConnectivity(), owner_set );
      if ( &edge_to_base != &edge_to_head ) {
         getNabrOwnerSet( edge_to_base.getConnectivity(), owner_set );
         getNabrOwnerSet( edge_to_base.d_partner->getConnectivity(), owner_set );
      }
      std::set<int>::iterator i_owner;
      for ( i_owner=owner_set.begin(); i_owner!=owner_set.end(); ++i_owner ) {
         if ( *i_owner != d_rank ) {
            mesg_to_owners[*i_owner];
         }
      }

   }

   const size_t number_of_mesg = mesg_to_owners.size();
   bridge_discoverEdges( edge_to_head,
                         edge_to_base,
                         mesg_to_owners,
                         d_partner == this );
   TBOX_ASSERT( number_of_mesg == mesg_to_owners.size() );
   TBOX_ASSERT( getNodeContainer().size() == getConnectivity().size() );
   TBOX_ASSERT( d_partner->getNodeContainer().size() == d_partner->getConnectivity().size() );
   if ( d_nprocs > 1 ) {
      bridge_shareEdges( mesg_to_owners );
      TBOX_ASSERT( getNodeContainer().size() == getConnectivity().size() );
      TBOX_ASSERT( d_partner->getNodeContainer().size() == d_partner->getConnectivity().size() );
   }

   return;
}



template<int DIM>
void LayerEdgeSet<DIM>::bridge_discoverEdges(
   const LayerEdgeSet &edge_to_head,
   const LayerEdgeSet &edge_to_base,
   std::map<LocalIndex, CommunicationStruct > &mesg_to_owners,
   const bool ignore_self_overlap )
{
   const NodeContainer &nodes = getNodeContainer();
   Connectivity &cnect = getConnectivity();


#if defined(DEBUG_SEND_EDGE_EXTRA_DATA) || defined(DEBUG_CHECK_ASSERTIONS)
   {
      /*
       * For assurance (debugging), put the message source, destination
       * and (yet unknown) size at the beginning of the messages.
       */
      typename std::map<LocalIndex,CommunicationStruct>::iterator i_mesg;
      for ( i_mesg=mesg_to_owners.begin();
            i_mesg!=mesg_to_owners.end();
            ++i_mesg ) {
         const int owner = (*i_mesg).first;
         CommunicationStruct &comm_struct = (*i_mesg).second;
         comm_struct.send_mesg.insert( comm_struct.send_mesg.end(),
                                       d_rank );
         comm_struct.send_mesg.insert( comm_struct.send_mesg.end(),
                                       owner );
         comm_struct.send_mesg.insert( comm_struct.send_mesg.end(),
                                       -1 );
      }
   }
#endif

   {
      /*
       * Ensure all neighbor lists exists.  In the loops to discover
       * edges, if an edge does not exist, the neighbor container
       * does not get created, ruining the one-to-one relationship
       * between nodes and neighbor container.
       *
       * This could be made much more efficient by replacing
       * dereferences with orderly insertion.
       */
      typename NodeContainer::const_iterator i_node;
      if ( d_parallel_state == GLOBALIZED ) {
         // d_cnect is the container for all processes.
         for ( i_node=nodes.begin(); i_node!=nodes.end(); ++i_node ) {
            d_cnect[(*i_node).getOwnerRank()][(*i_node).getLocalIndex()];
         }
      }
      else {
         for ( i_node=nodes.begin(); i_node!=nodes.end(); ++i_node ) {
            cnect[(*i_node).getLocalIndex()];
         }
      }
   }

   /*
    * Get from the two given legs the middle node layer
    * (common layer shared by the first and second legs),
    * the first layer (corresponding to internal layer)
    * and the second layer (corresponding to the target layer).
    */
   const NodeContainer &middle_nodes = edge_to_head.getNodeContainer();

   const hier::IntVector<DIM> &base_ratio = edge_to_base.getBaseRefinementRatio();
   const hier::IntVector<DIM> &head_ratio = edge_to_head.getBaseRefinementRatio();

   const bool base_is_finer = base_ratio > head_ratio;
   const hier::IntVector<DIM> base_head_ratio = base_is_finer ?
     base_ratio/head_ratio : hier::IntVector<DIM>(0);

   const bool head_is_finer = head_ratio > base_ratio;
   const hier::IntVector<DIM> head_base_ratio = head_is_finer ?
     head_ratio/base_ratio : hier::IntVector<DIM>(0);

   const int node_com_buffer_size = Node::commBufferSize();

   /*
    * Sanity checks.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( edge_to_base.getNodeContainer() == edge_to_head.getNodeContainer() );
#endif

   /*
     A neighbor of a node would be discarded if
     - the two are equal by comparison,
     - they are from layers with the same refinement ratio, and
     - ignore_self_overlap is true.
   */
   const bool discard_self_overlap =
      ignore_self_overlap &&
      ( d_base->getRefinementRatio() ==
        d_partner->d_base->getRefinementRatio() );

   /*
    * As the edges between the tail and head layers are discovered,
    * place the info directly into message buffers to be sent to
    * processes that need the info.  A dummy communication struct
    * is used in serial mode, where no message passing is needed.
    */
   std::vector<int> dummy_mesg;
   tbox::Array<int> tmp_buffer(2+Node::commBufferSize());

   /*
    * Use the middle layer, which has edges to the tail
    * and head layers of the bridge layer edge, to determine
    * which nodes to check for intersection, check them and
    * create edges from the tail to the head.  Send the edges
    * to the owners of the boxes in the tail.
    */
   typename NodeContainer::const_iterator i_midd;
   typename Connectivity::const_iterator i_nabr_to_base;
   typename Connectivity::const_iterator i_nabr_to_head;
   for ( i_midd=middle_nodes.begin(),
            i_nabr_to_base=edge_to_base.getConnectivity().begin(),
            i_nabr_to_head=edge_to_head.getConnectivity().begin();
         i_midd!=middle_nodes.end();
         ++i_midd,++i_nabr_to_base,++i_nabr_to_head ) {

      /*
       * If the head is at a finer or equal refinement ratio, 
       * grow the boxes in the head before checking intersection.
       * If the tail is finer, use the head as is, but grow each
       * box in the tail before checking.
       */
      // const NabrContainer &nabr_to_head = temp_nabr_to_head;
      const NabrContainer &nabr_to_head = (*i_nabr_to_head).second;
      const NabrContainer &nabr_to_base = (*i_nabr_to_base).second;

      int i;
      typename NabrContainer::const_iterator i_base;
      typename NabrContainer::const_iterator i_head;

      /*
       * Create a copy of the list of neighbors to the head,
       * so we can modify it for the purpose of comparison.
       * We either refine these nodes or grow them.  For
       * intersection checks, the boxes must be compared at
       * the same index space and the node on the finer layer
       * must be grown before checking.
       */
      hier::BoxArray<DIM> nabr_to_head_boxes((*i_nabr_to_head).second.size());
      if ( base_is_finer ) {
         for ( i=0,i_head=nabr_to_head.begin();
               i<nabr_to_head_boxes.size(); ++i,++i_head ) {
           nabr_to_head_boxes(i) = (*i_head).getBox();
           nabr_to_head_boxes(i).refine(base_head_ratio);
         }
      }
      else {
         for ( i=0,i_head=nabr_to_head.begin();
               i<nabr_to_head_boxes.size(); ++i,++i_head ) {
           nabr_to_head_boxes(i) = (*i_head).getBox();
           nabr_to_head_boxes(i).grow(d_partner->d_gcw);
         }
      }

      for ( i_base=nabr_to_base.begin();
            i_base!=nabr_to_base.end(); ++i_base ) {

         const Node &base_node = *i_base;
         std::vector<int> &mesg_to_base_owner =
            d_parallel_state == GLOBALIZED || base_node.getOwnerRank() == d_rank ?
            dummy_mesg : mesg_to_owners[base_node.getOwnerRank()].send_mesg;

         hier::Box<DIM> base_box = base_node.getBox();
         if ( base_is_finer ) base_box.grow(d_gcw);
         else if ( head_is_finer ) base_box.refine(head_base_ratio);

         /*
          * Loop through possible nodes at the head that may intersect.
          * Using hier::BoxTree<DIM> would be more efficient, but hier::BoxTree<DIM>
          * does not currently support Node objects yet
          * (just boxes).
          *
          * // hier::BoxList<DIM> overlap_boxes;
          * // head_boxtree.findOverlapBoxes( overlap_boxes, base_box );
          */
         for ( i=0,i_head=nabr_to_head.begin();
               i<nabr_to_head_boxes.size(); ++i,++i_head ) {

            const Node &head_node = *i_head;

            std::vector<int> &mesg_to_head_owner =
               d_parallel_state == GLOBALIZED || head_node.getOwnerRank() == d_rank ?
               dummy_mesg : mesg_to_owners[head_node.getOwnerRank()].send_mesg;
            /*
             * mesg_to_head_owner MUST be declared here instead of in
             * the following conditional block.  It creates the message,
             * which is important even if the message is empty.
             * We differentiate an empty message from no message
             * for coordinating communications.
             */

            if ( discard_self_overlap && head_node == base_node ) continue;

            if ( base_box.intersects(nabr_to_head_boxes(i)) ) {

               if ( d_parallel_state == GLOBALIZED ) {
                  d_cnect[base_node.getOwnerRank()]
                     [base_node.getLocalIndex()].insert(head_node);
                  if ( d_partner != this ) {
                     d_partner->d_cnect[head_node.getOwnerRank()]
                        [head_node.getLocalIndex()].insert(base_node);
                  }
               }
               else {

                  /*
                   * In the messages,
                   * each edge is represented by a string of integers.
                   * The first integer is 0 (representing information for
                   * the (forward) bridge from self to partner) or 1
                   * (representing information for the reverse bridge
                   * from partner to self).  The second integer is the
                   * local index of the node for whom the neighbor is found.
                   * The remaining integers are the neighbor data.
                   */

                  if ( base_node.getOwnerRank() == d_rank ) {
                     TBOX_ASSERT( cnect.find(base_node.getLocalIndex()) != cnect.end() );
                     cnect[base_node.getLocalIndex()].insert(head_node);
                  }
                  else {
                     tmp_buffer[0] = 0;
                     tmp_buffer[1] = base_node.getLocalIndex();
                     head_node.putToIntBuffer( &tmp_buffer[2] );
                     mesg_to_base_owner.insert( mesg_to_base_owner.end(),
                                                tmp_buffer.getPointer(),
                                                tmp_buffer.getPointer()+
                                                2+node_com_buffer_size );
                  }

                  if ( d_partner != this ) {
                     if ( head_node.getOwnerRank() == d_rank ) {
                        TBOX_ASSERT( d_partner->getConnectivity().find(head_node.getLocalIndex()) != d_partner->getConnectivity().end() );
                        d_partner->
                           getConnectivity()[head_node.getLocalIndex()].insert(base_node);
                     }
                     else {
                        tmp_buffer[0] = 1;
                        tmp_buffer[1] = head_node.getLocalIndex();
                        base_node.putToIntBuffer( &tmp_buffer[2] );
                        mesg_to_head_owner.insert( mesg_to_head_owner.end(),
                                                   tmp_buffer.getPointer(),
                                                   tmp_buffer.getPointer()+
                                                   2+node_com_buffer_size );
                     }
                  }
               }

            } // Intersection block.

         } // Small head-node loop.

      } // Small base-node loop.

   } // Middle-node loop.

#ifdef DEBUG_CHECK_ASSERTIONS
   checkNodeNabrCorrespondance();
#endif

#if defined(DEBUG_SEND_EDGE_EXTRA_DATA) || defined(DEBUG_CHECK_ASSERTIONS)
   {
      /*
       * For assurance (debugging), put the message
       * size near the beginning of the messages.
       */
      typename std::map<LocalIndex,CommunicationStruct>::iterator i_mesg;
      for ( i_mesg=mesg_to_owners.begin();
            i_mesg!=mesg_to_owners.end();
            ++i_mesg ) {
         CommunicationStruct &comm_struct = (*i_mesg).second;
         comm_struct.send_mesg[2] = comm_struct.send_mesg.size();
      }
   }
#endif

   return;
}



template<int DIM>
void LayerEdgeSet<DIM>::bridge_shareEdges(
  std:: map<LocalIndex,CommunicationStruct> &mesg_to_owners )
{

   // Nothing to do for serial run.
#if defined(HAVE_MPI)

   const int node_com_buffer_size = Node::commBufferSize();
   Connectivity &cnect = getConnectivity();

   /*
    * Send and receive the data-sharing message-lengths.
    * Allocate space for incoming data-sharing messages and post
    * send/receives for them.
    * While messages fly, unpack locally found edges that should
    * stay local.
    * Reveive incoming data-sharing messages, unpacking them as
    * they arrive.
    */
   typename std::map<LocalIndex,CommunicationStruct>::iterator i_mesg_to_owner;
   for ( i_mesg_to_owner=mesg_to_owners.begin();
         i_mesg_to_owner!=mesg_to_owners.end();
         ++i_mesg_to_owner ) {
      const int owner = (*i_mesg_to_owner).first;
      CommunicationStruct &comm_struct = (*i_mesg_to_owner).second;
      if ( owner != d_rank ) {
         comm_struct.recv_done = false;
         comm_struct.recv_size = -1;
         MPI_Irecv( &comm_struct.recv_size,
                    1,
                    MPI_INT,
                    owner,
                    SIZE_EXCHANGE_TAG,
                    tbox::SAMRAI_MPI::commWorld,
                    &comm_struct.recv_rqst );
         comm_struct.send_size = comm_struct.send_mesg.size();
         MPI_Isend( &comm_struct.send_size,
                    1,
                    MPI_INT,
                    owner,
                    SIZE_EXCHANGE_TAG,
                    tbox::SAMRAI_MPI::commWorld,
                    &comm_struct.send_rqst );
         MPI_Request_free( &comm_struct.send_rqst );
         tbox::plog << "shareEdges proc " << d_rank << " sends " << comm_struct.send_size << " to proc " << owner << std::endl;
         if ( comm_struct.send_size > 0 ) {
            comm_struct.send_done = false;
            MPI_Isend( &comm_struct.send_mesg[0],
                       comm_struct.send_mesg.size(),
                       MPI_INT,
                       owner,
                       EDGE_EXCHANGE_TAG,
                       tbox::SAMRAI_MPI::commWorld,
                       &comm_struct.send_rqst );
         }
         else {
            comm_struct.send_done = true;
         }
      }
      else {
         comm_struct.send_done = comm_struct.recv_done = true;
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( mesg_to_owners.find(d_rank) == mesg_to_owners.end() );
#endif
   int n_outstanding_mesg = mesg_to_owners.size();
   i_mesg_to_owner=mesg_to_owners.begin();
   while ( n_outstanding_mesg > 0 ) {
      const int owner = (*i_mesg_to_owner).first;
      CommunicationStruct &comm_struct = (*i_mesg_to_owner).second;
      ++i_mesg_to_owner;
      if ( i_mesg_to_owner == mesg_to_owners.end() ) {
         i_mesg_to_owner = mesg_to_owners.begin();
      }
      if ( !comm_struct.recv_done ) {
         int flag;
         tbox::SAMRAI_MPI::status status;
         MPI_Test( &comm_struct.recv_rqst, &flag, &status );
         if ( flag && status.MPI_TAG == SIZE_EXCHANGE_TAG ) {
            // tbox::plog << "shareEdges proc " << d_rank << " expects " << comm_struct.recv_size << " from proc " << owner << endl;
            int tmp_count = -1;
            MPI_Get_count( &status, MPI_INT, &tmp_count );
            TBOX_ASSERT( tmp_count == 1 );
            if ( comm_struct.recv_size == 0 ) {
               --n_outstanding_mesg;
               comm_struct.recv_done = true;
            }
            else {
               /*
                * Allocate space and post receive for edge info as soon
                * as we know the size of the incoming message.
                */
               comm_struct.recv_mesg.clear();
               comm_struct.recv_mesg.insert( comm_struct.recv_mesg.end(),
                                             comm_struct.recv_size,
                                             0 );
               MPI_Irecv( &comm_struct.recv_mesg[0],
                          comm_struct.recv_size,
                          MPI_INT,
                          owner,
                          EDGE_EXCHANGE_TAG,
                          tbox::SAMRAI_MPI::commWorld,
                          &comm_struct.recv_rqst );
            }
         }
         else if ( flag && status.MPI_TAG == EDGE_EXCHANGE_TAG ) {
            TBOX_ASSERT( comm_struct.recv_mesg.size() > 0 ); // Means size already received.
            int tmp_count = -1;
            MPI_Get_count( &status, MPI_INT, &tmp_count );
            TBOX_ASSERT( status.MPI_SOURCE == owner );
            TBOX_ASSERT( status.MPI_TAG == EDGE_EXCHANGE_TAG );
            TBOX_ASSERT( tmp_count == comm_struct.recv_size );
            /*
             * Unpack message as soon as it arrives and set the done flag.
             */
            --n_outstanding_mesg;
            comm_struct.recv_done = true;
            const std::vector<int> &recv_mesg = comm_struct.recv_mesg;
            std::vector<int>::const_iterator i_recv_mesg = recv_mesg.begin();
#if defined(DEBUG_SEND_EDGE_EXTRA_DATA) || defined(DEBUG_CHECK_ASSERTIONS)
            {
               /*
                * For assurance (debugging), get the message source,
                * destination and size at the beginning of the messages.
                */
               int tmp_src = *(i_recv_mesg++);
               int tmp_dst = *(i_recv_mesg++);
               int tmp_sze = *(i_recv_mesg++);
               TBOX_ASSERT( tmp_src == owner );
               TBOX_ASSERT( tmp_dst == d_rank );
               TBOX_ASSERT( tmp_sze == comm_struct.recv_size );
            }
#endif
            tbox::plog << "shareEdges proc " << d_rank << " recvs " << comm_struct.recv_size << " from proc " << owner << std::endl;
            for ( ;
                  i_recv_mesg!=recv_mesg.end();
                  i_recv_mesg+=(2+node_com_buffer_size) ) {
               const int *buffer = &(*i_recv_mesg);
               const int h_or_t = buffer[0];
               const int local_node_index = buffer[1];
               Node nabr;
               nabr.getFromIntBuffer(buffer+2);
               if ( h_or_t == 0 ) {
                  // This neighbor is node local_node_index in the tail layer.
                  TBOX_ASSERT( cnect.find(local_node_index) != cnect.end() );
                  cnect[local_node_index].insert(
                    cnect[local_node_index].end(), nabr );
               }
               else {
                  // This neighbor is node local_node_index in the head layer.
                  TBOX_ASSERT( d_partner->getConnectivity().find(local_node_index) != d_partner->getConnectivity().end() );
                  d_partner->getConnectivity()[local_node_index].insert(
                    d_partner->getConnectivity()[local_node_index].end(), nabr );
               }
            }
         }
      }
   }

   /*
    * Free the send requests and as assurance,
    * check that they are completed without errors.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( mesg_to_owners.find(d_rank) == mesg_to_owners.end() );
#endif
   n_outstanding_mesg = mesg_to_owners.size();
   i_mesg_to_owner=mesg_to_owners.begin();
   while ( n_outstanding_mesg > 0 ) {
      CommunicationStruct &comm_struct = (*i_mesg_to_owner).second;
      ++i_mesg_to_owner;
      if ( i_mesg_to_owner == mesg_to_owners.end() ) {
         i_mesg_to_owner = mesg_to_owners.begin();
      }
      if ( !comm_struct.send_done ) {
         int flag;
         tbox::SAMRAI_MPI::status status;
         MPI_Test( &comm_struct.send_rqst, &flag, &status );
         if ( flag ) {
            --n_outstanding_mesg;
            comm_struct.send_done = true;
         }
      }
   }

#endif

   return;
}






/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::setParallelState( const ParallelState parallel_state )
{

#ifdef REQUIRE_PARTNER_PARALLEL_STATE_MATCH
   bool work_on_partner = ( d_partner && d_partner != this );
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( work_on_partner ) {
     TBOX_ASSERT ( d_partner->d_parallel_state == d_parallel_state );
   }
#endif
#endif

   int n;
   Connectivity *tmp_cnect = NULL;

   if ( parallel_state != d_parallel_state ) {


      if ( d_parallel_state == DISTRIBUTED &&
           parallel_state == GLOBALIZED ) {
         // Going to GLOBALIZED state.
         tmp_cnect = new Connectivity[d_nprocs];
         TBOX_ASSERT( tmp_cnect != NULL );
         tmp_cnect[d_rank] = *d_cnect;
         d_cnect->clear();
         d_cnect = tmp_cnect;
         tmp_cnect = NULL;
#ifdef REQUIRE_PARTNER_PARALLEL_STATE_MATCH
         if ( work_on_partner ) {
            tmp_cnect = new Connectivity[d_nprocs];
            TBOX_ASSERT( tmp_cnect != NULL );
            d_partner->d_cnect[d_rank] = *d_partner->d_cnect;
            d_partner->d_cnect->clear();
            d_partner->d_cnect = tmp_cnect;
            tmp_cnect = NULL;
         }
#endif
         for ( n=0; n<d_nprocs; ++n ) if ( n != d_rank ) {
           allocateConnectivity( getNodeContainer(n), d_cnect[n] );
#ifdef REQUIRE_PARTNER_PARALLEL_STATE_MATCH
           if ( work_on_partner ) {
              allocateConnectivity( d_partner->getNodeContainer(n),
                                    d_partner->d_cnect[n] );
           }
#endif
         }
      }


      else if ( d_parallel_state == GLOBALIZED &&
                parallel_state == DISTRIBUTED ) {
         // Going to DISTRIBUTED state.
         if ( d_base->getParallelState() != LayerNodeSet<DIM>::DISTRIBUTED ) {
            TBOX_ERROR("Attempt to set GLOBALIZED mode on a LayerEdgeSet\n"
                       <<"without first setting GLOBALIZED mode on its\n"
                       <<"LayerNodeSet.  The LayerEdgeSet cannot be\n"
                       <<"GLOBALIZED when its LayerNodeSet is not.\n"
                       <<"The LayerNodeSet is an external object which\n"
                       <<"which cannot be modified by the LayerEdgeSet.\n");
         }
         tmp_cnect = new Connectivity[1];
         TBOX_ASSERT( tmp_cnect != NULL );
         *tmp_cnect = d_cnect[d_rank];
         delete [] d_cnect;
         d_cnect = tmp_cnect;
         tmp_cnect = NULL;
#ifdef REQUIRE_PARTNER_PARALLEL_STATE_MATCH
         if ( work_on_partner ) {
            tmp_cnect = new Connectivity[1];
            TBOX_ASSERT( tmp_cnect != NULL );
            *tmp_cnect = d_partner->d_cnect[d_rank];
            delete [] d_partner->d_cnect;
            d_partner->d_cnect = tmp_cnect;
            tmp_cnect = NULL;
         }
#endif
      }

      d_parallel_state = parallel_state;
#ifdef REQUIRE_PARTNER_PARALLEL_STATE_MATCH
      if ( work_on_partner ) {
        d_partner->d_parallel_state = parallel_state;
      }
#endif

   }
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::initialize(
   ParallelState parallel_state,
   const LayerNodeSet<DIM> &layer_node_set,
   const hier::IntVector<DIM> &head_refinement_ratio,
   const hier::IntVector<DIM> &gcw,
   const Connectivity *connectivity )
{
   if ( d_cnect != NULL ) {
      delete [] d_cnect;
   }

   if ( !(head_refinement_ratio > hier::IntVector<DIM>(0)) ) {
      TBOX_ERROR("LayerEdgeSet::initialize():\n"
                 <<"Invalid refinement ratio given: "
                 << head_refinement_ratio << "\n");
   }
   d_head_ratio = head_refinement_ratio;

   if ( gcw < hier::IntVector<DIM>(0) ) {
      TBOX_ERROR("LayerEdgeSet::initialize():\n"
                 <<"Invalid ghost cell width: "
                 << gcw << "\n");
   }
   d_gcw = gcw;

   int n, nbeg=-1, nend=-1;

   switch (parallel_state) {
   case DISTRIBUTED:
      nbeg = d_rank;
      nend = d_rank+1;
      d_cnect = new Connectivity [1];
      break;
   case GLOBALIZED:
      nbeg = 0;
      nend = d_nprocs;
      if ( layer_node_set.getParallelState() != LayerNodeSet<DIM>::GLOBALIZED ) {
         TBOX_ERROR("Setting GLOBALIZED connectivity requires the layer node set\n"
                    "to be in GLOBALIZED mode first.\n");
      }
      d_cnect = new Connectivity [d_nprocs];
      break;
   default:
      TBOX_ERROR("Unrecognized parallel_state '" << parallel_state << "'\n");
   }

   if ( connectivity ) {

      for ( n=nbeg; n<nend; ++n ) {

         const Connectivity &cnect =
            connectivity[parallel_state==DISTRIBUTED?0:n];

#ifdef DEBUG_CHECK_ASSERTIONS
         const NodeContainer &node_container =
            layer_node_set.getNodeContainer(n);
         if ( node_container.size() != cnect.size() ) {
            TBOX_ERROR("Given node layer set and connectivity have\n"
                       <<"incompatible size for process " << n << "\n");
         }
         typename NodeContainer::const_iterator ni;
         typename Connectivity::const_iterator ci;
         for ( ni=node_container.begin(); ni!=node_container.end(); ++ni ) {
            ci = cnect.find( (*ni).getLocalIndex() );
            if ( ci == cnect.end() ) {
               TBOX_ERROR("Node " << (*ni) << " appears in the layer node set\n"
                          <<"but not in the given connecvity.\n");
            }
         }
#endif

         d_cnect[n] = cnect;

      }
   }

   d_parallel_state = parallel_state;
   d_base = &layer_node_set;

   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::getNabrOwnerSet( const Connectivity &cnect,
                                           std::set<int> &owners ) const
{
   typename Connectivity::const_iterator i_node;
   for ( i_node=cnect.begin(); i_node!=cnect.end(); ++i_node ) {
      const typename Connectivity::value_type &node = *i_node;
      const NabrContainer &nabr = node.second;
      typename NabrContainer::const_iterator i_nabr;
      for ( i_nabr=nabr.begin(); i_nabr!=nabr.end(); ++i_nabr ) {
         const int owner = (*i_nabr).getOwnerRank();
         owners.insert(owner);
      }
   }
   return;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::printClassData( std::ostream &co, int detail_depth ) const
{
   const Connectivity &cnect = getConnectivity();
   const NodeContainer &nodes = getNodeContainer();
   co << "Parallel state   : " << (getParallelState()==DISTRIBUTED?"DIST":"GLOB") << '\n'
      << "Base,head ratio  : " << getBaseRefinementRatio() << ", " << getHeadRefinementRatio() << '\n'
      << "Ghost cell width : " << d_gcw << '\n'
      << "Node count       : " << nodes.size() << ", " << cnect.size() << '\n'
      ;
   if ( detail_depth > 0 ) {
      typename NodeContainer::const_iterator i_node;
      typename Connectivity::const_iterator i_node_nabr;
      co << "Nodes:\n";
      for ( i_node=nodes.begin(), i_node_nabr=cnect.begin();
            i_node!=nodes.end();
            ++i_node ) {
         co << "    "
            << (*i_node) << "\t"
            << (*i_node).getBox().numberCells() << '\n';
         if ( i_node_nabr != cnect.end() ) {
            TBOX_ASSERT( (*i_node).getLocalIndex() == (*i_node_nabr).first );
            const NabrContainer &nabrs = (*i_node_nabr).second;
            co << "    Nabrs (" << nabrs.size() << "):\n";
            if ( detail_depth > 1 ) {
               typename NabrContainer::const_iterator i_nabr;
               for ( i_nabr=nabrs.begin(); i_nabr!=nabrs.end(); ++i_nabr ) {
                  co << "        "
                     << (*i_nabr) << "\t"
                     << (*i_nabr).getBox().numberCells() << '\n';
               }
            }
            else {
               co << " ...\n";
            }
            ++i_node_nabr;
         }
      }
   }
   return;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::printEdgeStats( std::ostream &co ) const
{
   const NodeContainer &nodes = getNodeContainer();
   const Connectivity &cnect = getConnectivity();
   /*
    * Compute additional statistics.
    * i_node iterates through the node list.
    * i_nabr iterates through the node-neighbor list,
    *        which is either a single list in non-serial mode
    *        or a series of lists in serial mode.  Of the seies
    *        of lists, node_nabr points to the current list being
    *        iterated through.
    */
   typename Connectivity::const_iterator i_nabr = d_parallel_state == GLOBALIZED ?
      d_cnect->begin() : cnect.begin();
   typename NodeContainer::const_iterator i_node;
   Connectivity *node_nabr = d_cnect;
   int tot_nnabr=0, min_nnabr=999999, max_nnabr=0;
   double avg_nnabr;
   for ( i_node=nodes.begin(); i_node!=nodes.end(); ++i_node ) {
      const int nnabr = (*i_nabr).second.size();
      tot_nnabr += nnabr;
      if ( nnabr < min_nnabr ) min_nnabr = nnabr;
      if ( nnabr > max_nnabr ) max_nnabr = nnabr;
      ++i_nabr;
      if ( d_parallel_state == GLOBALIZED && i_nabr == node_nabr->end() &&
           node_nabr < d_cnect+d_nprocs-1 ) {
         ++node_nabr;
         i_nabr = node_nabr->begin();
      }
   }
   avg_nnabr = tot_nnabr == 0 ? 0.0 : double(tot_nnabr)/nodes.size();

   co << "N nabrs (l,h,t,a): "
      << std::setw(10) << min_nnabr
      << std::setw(10) << max_nnabr
      << std::setw(10) << tot_nnabr
      << std::setw(10) << avg_nnabr << '\n'
      ;

   std::set<int> semilocal_owners;
   getNabrOwnerSet( getConnectivity(), semilocal_owners );
   co << "N of semilocal owners: " << semilocal_owners.size() << '\n';
   
   return;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
LayerNodeSet<DIM> *LayerEdgeSet<DIM>::copyAndGlobalize( const LayerNodeSet<DIM> &r )
{
   // Prevent wasteful accidental use when this method is not needed.
   TBOX_ASSERT( r.getParallelState() != LayerNodeSet<DIM>::GLOBALIZED );

   LayerNodeSet<DIM> *copy = new LayerNodeSet<DIM>(r);

   copy->setParallelState(LayerNodeSet<DIM>::GLOBALIZED);

   return copy;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::checkNodeNabrCorrespondance() const {
   const NodeContainer &nodes = getNodeContainer();
   const Connectivity &cnect = getConnectivity();
   /*
    * There should a one-to-one correspondance between nodes and
    * neighbor containers.  It is mainly for debugging this class.
    */
   if ( d_parallel_state == GLOBALIZED ) {
      // Checking can be done here, too, but I've not written it.
      TBOX_WARNING("checkNodeNabrCorrespondance() not yet written\n"
                   << "for GLOBALIZED mode.\n");
   }
   else {
      if ( nodes.size() != cnect.size() ) {
         TBOX_ERROR("Mismatched nodes and neighbor container sizes.\n"
                    << nodes.size() << " nodes vs "
                    << cnect.size() << " neighbor containers.\n"
                    << "This is an error in the library!");
      }
      typename NodeContainer::const_iterator i_node;
      typename Connectivity::const_iterator i_nabr;
      for ( i_node=nodes.begin(), i_nabr=cnect.begin();
            i_node!=nodes.end(); ++i_node, ++i_nabr ) {
         if ( (*i_node).getLocalIndex() != (*i_nabr).first ) {
            TBOX_ERROR("Mismatched nodes and neighbor container.\n"
                       << "Node " << (*i_node).getLocalIndex()
                       << "lines up with container " << (*i_nabr).first << "\n"
                       << "This is an error in the library!");
         }
      }
   }
   return;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::checkConnectivity( const LayerNodeSet<DIM> &head ) const
{

   /*
    * Rebuild the edge containers, then check that
    * the rebuilt edges match the existing edges.
    *
    * Currently, we use findEdges to rebuild the edges.
    * Thus, it may be pointless to use this method
    * as a check for that method.
    */
   LayerEdgeSet rebuilt;
   rebuilt.initialize( DISTRIBUTED, *d_base, d_gcw );

   /*
    * To rebuild the layer without help, we need the global boxes
    * for the layer edge's head (the partner's base).
    * Use head if it is GLOBALIZED;
    * otherwise, build a temporary one with the global boxes.
    * (We are not allowed to modify d_partner->d_base
    * because it is managed externally!)
    */
   if ( head.getParallelState() ==
        LayerNodeSet<DIM>::GLOBALIZED ) {
      rebuilt.findEdges_rbbt( head, d_base == &head );
   }
   else {
      LayerNodeSet<DIM> *temp_head =
         copyAndGlobalize(head);
      rebuilt.findEdges_rbbt( *temp_head, d_base == &head );
      delete temp_head;
   }

   {
      /*
       * Check the current data against the rebuilt data.
       */

      const Connectivity &tcnect = *d_cnect;
      const Connectivity &rcnect = *rebuilt.d_cnect;
      if ( tcnect.size() != rcnect.size() ) {
         TBOX_ERROR("Wrong number of neighbor lists\n");
      }

      const NodeContainer &tnodes = getNodeContainer();
      const NodeContainer &rnodes = rebuilt.getNodeContainer();
      typename NodeContainer::const_iterator i_tnode, i_rnode;

      for ( i_tnode=tnodes.begin(), i_rnode=rnodes.begin();
            i_tnode!=tnodes.end();
            ++i_tnode, ++i_rnode ) {

         const Node &tnode = *i_tnode;
         const Node &rnode = *i_rnode;
         typename Connectivity::const_iterator
            i_tnabrs = tcnect.find( tnode.getLocalIndex() ),
            i_rnabrs = rcnect.find( rnode.getLocalIndex() );

         const NabrContainer &tnabrs = (*i_tnabrs).second;
         const NabrContainer &rnabrs = (*i_rnabrs).second;
         if ( tnabrs != rnabrs ) {
            typename NabrContainer::const_iterator ni;
            tbox::pout << "Neighbors for " << tnode << " (" << rnode << "):\n";
            tbox::pout << "true:\n";
            for ( ni=tnabrs.begin(); ni!=tnabrs.end(); ++ni ) {
               tbox::pout << *ni << std::endl;
            }
            tbox::pout << "rebuilt:\n";
            for ( ni=rnabrs.begin(); ni!=rnabrs.end(); ++ni ) {
               tbox::pout << *ni << std::endl;
            }
            TBOX_ERROR("Wrong edge data for node "
                       << (*i_tnode).getLocalIndex() << "\n");
         }
      }
   }

   return;
}




/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::checkNodeConsistency() const {
   if ( d_partner == NULL ) {
      TBOX_ERROR("Attached partner is required in checkNodeConsistency().\n");
   }

   LayerNodeSet<DIM> *temp_head = NULL;
   if ( d_partner->d_base->getParallelState() ==
        LayerNodeSet<DIM>::DISTRIBUTED ) {
      temp_head = copyAndGlobalize(*d_partner->d_base);
   }

   checkNodeConsistency( getNodeContainer(),
                         getConnectivity(),
                         *( temp_head ? (const LayerNodeSet<DIM>*)temp_head :
                            d_partner->d_base ) );

   if ( temp_head != NULL ) {
      delete temp_head;
   }

   return;
}

/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerEdgeSet<DIM>::checkNodeConsistency(
   const NodeContainer &base_node_container,
   const Connectivity &cnect,
   const LayerNodeSet<DIM> &head_layer_node_set )
{
   
   typename NodeContainer::const_iterator i_local_nodes;
   for ( i_local_nodes=base_node_container.begin();
         i_local_nodes!=base_node_container.end();
         ++i_local_nodes ) {
      /*
        For each node in the given base_node_container,
        check that there is a neighbor list for it.
        For each neighbor in each neighbor list,
        check that the neighbor is in the head_layer_node_set.
      */

      const Node &node = *i_local_nodes;

      typename Connectivity::const_iterator i_cnect =
         cnect.find(node.getLocalIndex());
      if ( i_cnect == cnect.end() ) {
         TBOX_ERROR("Missing connectivity data for node "
                    << node << "\n");
      }

      const NabrContainer &nabrs = (*i_cnect).second;

      typename NabrContainer::const_iterator i_nabr;
      for ( i_nabr=nabrs.begin(); i_nabr!=nabrs.end(); ++i_nabr ) {

         const Node &nabr = *i_nabr;

         const NodeContainer &head_nodes =
            head_layer_node_set.getNodeContainer(nabr.getOwnerRank());

         typename NodeContainer::const_iterator i_nabr_in_head =
            head_nodes.find(nabr);

         if ( i_nabr_in_head == head_nodes.end() ) {
            TBOX_ERROR("Missing connectivity data at node "
                       << node << "\n");
         }

         const Node &nabr_in_head = *i_nabr_in_head;
         if ( nabr != nabr_in_head ||
              nabr.getBox() != nabr_in_head.getBox() ) {
            TBOX_ERROR("Inconsistent node data at node "
                       << node << "\n"
                       <<"Neighbor " << nabr << "does not match\n"
                       <<"node " << nabr_in_head << "\n");
         }

      }
   }

   return;
}


}
}
#endif
