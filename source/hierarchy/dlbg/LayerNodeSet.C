/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/dlbg/LayerNodeSet.C $
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1846 $
 * Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
 * Description: Set of layer nodes in a distributed box graph.
 */

#ifndef included_hier_LayerNodeSet_C
#define included_hier_LayerNodeSet_C

#include "LayerNodeSet.h"

#include "tbox/SAMRAI_MPI.h"

#include IOSTREAM_HEADER_FILE
#include IOMANIP_HEADER_FILE


#ifdef DEBUG_NO_INLINE
#include "LayerNodeSet.I"
#endif

namespace SAMRAI {
namespace hier {


#define SIZE_EXCHANGE_TAG (100000)
#define EDGE_EXCHANGE_TAG (100001)


template<int DIM>
LayerNodeSet<DIM>::LayerNodeSet()
   : d_nodes(NULL),
     d_ratio(0),
     d_parallel_state(DISTRIBUTED),
     d_rank(tbox::SAMRAI_MPI::getRank())
{
   d_nodes = new NodeContainer[1];


   t_acquire_nonlocal_nodes = tbox::TimerManager::getManager()->
      getTimer("LayerNodeSet::acquireNonlocalNodes()");

   TBOX_ASSERT( d_nodes != NULL );

   return;
}


template<int DIM>
LayerNodeSet<DIM>::LayerNodeSet( const LayerNodeSet &r)
   : d_nodes(NULL),
     d_ratio(r.d_ratio),
     d_parallel_state(r.d_parallel_state),
     d_rank(r.d_rank)
{
   if ( d_parallel_state == DISTRIBUTED ) {
      d_nodes = new NodeContainer[1];
      *d_nodes = *r.d_nodes;
   }
   else if ( d_parallel_state == GLOBALIZED ) {
      int nprocs = tbox::SAMRAI_MPI::getNodes();
      d_nodes = new NodeContainer[tbox::SAMRAI_MPI::getNodes()];
      int n;
      for ( n=0; n<nprocs; ++n ) {
         d_nodes[n] = r.d_nodes[n];
      }
   }


   t_acquire_nonlocal_nodes = tbox::TimerManager::getManager()->
      getTimer("LayerNodeSet::acquireNonlocalNodes()");


   TBOX_ASSERT( d_nodes != NULL );

   return;
}


template<int DIM>
LayerNodeSet<DIM>::~LayerNodeSet()
{
   deallocateData();
   delete [] d_nodes;
   d_nodes = NULL;
   return;
}



template<int DIM>
LayerNodeSet<DIM> &LayerNodeSet<DIM>::operator=( const LayerNodeSet &r )
{
   d_parallel_state = r.d_parallel_state;
   d_ratio = r.d_ratio;
   if ( d_parallel_state == DISTRIBUTED ) {
      d_nodes = new NodeContainer[1];
      *d_nodes = *r.d_nodes;
   }
   else {
      const int nprocs = tbox::SAMRAI_MPI::getNodes();
      d_nodes = new NodeContainer[nprocs];
      int n;
      for ( n=0; n<nprocs; ++n ) {
         d_nodes[n] = r.d_nodes[n];
      }
   }
   return *this;
}


template<int DIM>
const typename LayerNodeSet<DIM>::NodeContainer
&LayerNodeSet<DIM>::getNodeContainer( const int rank ) const
{
   if ( d_parallel_state == DISTRIBUTED && rank > -1 && rank != d_rank ) {
      TBOX_ERROR("Non-local nodes are not available in DISTRIBUTED distribution mode.");
   }
   return d_parallel_state == GLOBALIZED ? d_nodes[rank>-1?rank:d_rank] : *d_nodes;
}


template<int DIM>
typename LayerNodeSet<DIM>::NodeContainer
&LayerNodeSet<DIM>::getMutableNodeContainer( const int rank )
{
   if ( d_parallel_state == DISTRIBUTED && rank > -1 && rank != d_rank ) {
      TBOX_ERROR("Non-local nodes are not available in DISTRIBUTED distribution mode.");
   }
   return d_parallel_state == GLOBALIZED ? d_nodes[rank>-1?rank:d_rank] : *d_nodes;
}


template<int DIM>
const hier::IntVector<DIM> &LayerNodeSet<DIM>::getRefinementRatio() const
{
   return d_ratio;
}


template<int DIM>
void LayerNodeSet<DIM>::setRefinementRatio( const hier::IntVector<DIM> &ratio )
{
   d_ratio = ratio;
   return;
}


template<int DIM>
void LayerNodeSet<DIM>::deallocateData()
{
   const int nproc = tbox::SAMRAI_MPI::getNodes();
   int n;
   switch (d_parallel_state) {
   case DISTRIBUTED:
      d_nodes->clear();
      break;
   case GLOBALIZED:
      for ( n=0; n<nproc; ++n ) {
         d_nodes[n].clear();
      }
      break;
   default:
      TBOX_ERROR("Corrupted distribution state.\n");
   }
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::setParallelState( const ParallelState parallel_state )
{
   if ( parallel_state != DISTRIBUTED && parallel_state != GLOBALIZED ) {
      TBOX_ERROR("LayerNodeSet::setParallelState: Invalid distribution state: "
                 << parallel_state << "\n");
   }

   const int nproc = tbox::SAMRAI_MPI::getNodes();
   NodeContainer *tmp_nodes = NULL;

   if ( parallel_state != d_parallel_state ) {


      if ( d_parallel_state == DISTRIBUTED &&
           parallel_state == GLOBALIZED ) {
         // Going to more GLOBALIZED state.
         tmp_nodes = new NodeContainer[nproc];
         TBOX_ASSERT( tmp_nodes != NULL );
         tmp_nodes[d_rank] = *d_nodes;
         delete [] d_nodes;
         d_nodes = tmp_nodes;
         tmp_nodes = NULL;
         acquireNonlocalNodes();
      }


      else if ( d_parallel_state == GLOBALIZED &&
                parallel_state == DISTRIBUTED ) {
         // Going to DISTRIBUTED state.
         tmp_nodes = new NodeContainer[1];
         TBOX_ASSERT( tmp_nodes != NULL );
         *tmp_nodes = d_nodes[d_rank];
         delete [] d_nodes;
         d_nodes = tmp_nodes;
         tmp_nodes = NULL;
      }

      d_parallel_state = parallel_state;

   }
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::setParallelState( LayerNodeSet *layers,
                                      const int num_layers,
                                      const ParallelState parallel_state ) const
{
   if ( parallel_state != DISTRIBUTED && parallel_state != GLOBALIZED ) {
      TBOX_ERROR("LayerNodeSet::setParallelState: Invalid distribution state: "
                 << parallel_state << "\n");
   }
   int i;
   for ( i=0; i<num_layers; ++i ) {
      layers[i].setParallelState(parallel_state);
   }
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::acquireNonlocalNodes( const int hogging_process )
{
   LayerNodeSet *object = this;
   acquireNonlocalNodes( 1, &object, hogging_process );
   return;
}





/*
***********************************************************************
* Acquire nonlocal nodes for multiple layer node sets.                *
* This method combines communication for the multiple layers          *
* to increase message passing efficiency.                             *
*                                                                     *
* Note: This method is stateless (could be static).                   *
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::acquireNonlocalNodes(
   const int num_sets,
   LayerNodeSet *layer_node_sets[],
   const int hogging_process )
{
   // Nothing to do for a serial run.
#if defined(HAVE_MPI)

   t_acquire_nonlocal_nodes->start();
   int n;
#ifdef DEBUG_CHECK_ASSERTIONS
   for ( n=0; n<num_sets; ++n ) {
      if ( layer_node_sets[n]->getParallelState() != DISTRIBUTED ) {
         TBOX_ERROR("LayerNodeSet objects must be in distributed mode\n"
                    <<"when acquiring nonlocal nodes.\n");
      }
      if ( hogging_process >= 0 &&
           d_rank !=  hogging_process &&
           layer_node_sets[n]->d_nodes[0].size() != 0 ) {
         TBOX_ERROR("Bad assumption of a hogging process.\n"
                    <<"You assumed that process " << hogging_process << '\n'
                    <<"is the only one with any nodes.  However, process\n"
                    <<d_rank << " has " << layer_node_sets[n]->d_nodes[0].size()
                    <<" nodes in the " << n << "th LayerNodeSet.\n"
                    <<"To continue would lead to lost or messages.\n");
      }
   }
#endif

   const int nproc = tbox::SAMRAI_MPI::getNodes();

   /*
    * The presence of a hogging process indicates that no other
    * process has any data to share.
    * If there is a hogging process, we broadcast from there.
    * If not, do an all-gather.
    */
   if ( hogging_process < 0 ) {

      tbox::Array<int> send_mesg;
      tbox::Array<int> recv_mesg;
      /*
       * Pack nodes from all layer node sets into a single message.
       * Note that each layer node set object packs the size of its
       * sub-message into send_mesg.
       */
      for ( n=0; n<num_sets; ++n ) {
         const LayerNodeSet &layer_node_set = *layer_node_sets[n];
         layer_node_set.acquireNonlocalNodes_pack(send_mesg,
                                                  send_mesg.getSize());
      }
      int send_mesg_size = send_mesg.getSize();

      /*
       * Send and receive the data.
       */

      tbox::Array<int> recv_mesg_size(nproc);
      MPI_Allgather( &send_mesg_size,
                     1,
                     MPI_INT,
                     recv_mesg_size.getPointer(),
                     1,
                     MPI_INT,
                     tbox::SAMRAI_MPI::getCommunicator());

      tbox::Array<int> proc_offset(nproc);
      int totl_size = 0;
      for ( n=0; n<nproc; ++n ) {
         proc_offset[n] = totl_size;
         totl_size += recv_mesg_size[n];
      }
      recv_mesg.resizeArray(totl_size);
      MPI_Allgatherv( send_mesg.getPointer(),
                      send_mesg_size,
                      MPI_INT,
                      recv_mesg.getPointer(),
                      recv_mesg_size.getPointer(),
                      proc_offset.getPointer(),
                      MPI_INT,
                      tbox::SAMRAI_MPI::getCommunicator() );

      /*
       * Extract node info received from other processors.
       */
      for ( n=0; n<num_sets; ++n ) {
         LayerNodeSet &layer_node_set = *layer_node_sets[n];
         layer_node_set.acquireNonlocalNodes_unpack(recv_mesg, proc_offset);
      }

   }
   else {
      int mesg_size = -1;
      tbox::Array<int> mesg;
      if ( d_rank == hogging_process ) {
         for ( n=0; n<num_sets; ++n ) {
            const LayerNodeSet &layer_node_set = *layer_node_sets[n];
            layer_node_set.acquireNonlocalNodes_pack(mesg, mesg.getSize());
         }
         mesg_size = mesg.getSize();
      }
#if 1
      int one = 1;
      tbox::SAMRAI_MPI::bcast( &mesg_size, one, hogging_process );
#else
      MPI_Bcast( &mesg_size,
                 1,
                 MPI_INT,
                 hogging_process,
                 tbox::SAMRAI_MPI::getCommunicator());
#endif
      if ( d_rank != hogging_process ) {
         mesg.resizeArray(mesg_size);
      }
#if 1
      tbox::SAMRAI_MPI::bcast( mesg.getPointer(), mesg_size, hogging_process );
#else
      MPI_Bcast( mesg.getPointer(),
                 mesg_size,
                 MPI_INT,
                 hogging_process,
                 tbox::SAMRAI_MPI::getCommunicator());
#endif
      if ( d_rank != hogging_process ) {
         const int *ptr = mesg.getPointer();
         for ( n=0; n<num_sets; ++n ) {
            LayerNodeSet &layer_node_set = *layer_node_sets[n];
            ptr = layer_node_set.acquireNonlocalNodes_unpack(ptr, hogging_process);
         }
      }
   }

   t_acquire_nonlocal_nodes->stop();
#endif

   return;
}


/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::acquireNonlocalNodes_pack(
   tbox::Array<int> &send_mesg,
   int offset ) const
{
   NULL_USE(offset);
   const int rank = tbox::SAMRAI_MPI::getRank();
   /*
    * Node acquisition is done during moves away from distributed state.
    * Thus, do not rely on current value of d_parallel_state.  Instead note that
    * the meaning of d_nodes has been altered in setParallelState() but the
    * value of d_parallel_state has not been changed.
    *
    * Thus, d_nodes should already be allocated for hybrid mode.
    */
   const NodeContainer &nodes = d_nodes[rank];

   /*
    * Pack node info from d_nodes into send_mesg,
    * starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of nodes from self
    *   - Self nodes
    */
   const int node_com_buf_size = LayerNode<DIM>::commBufferSize();
   const int send_mesg_size = 1 + node_com_buf_size * ( nodes.size() );
   const int old_size = send_mesg.getSize();
   send_mesg.resizeArray(old_size + send_mesg_size);

   int *ptr = send_mesg.getPointer() + old_size;
   *(ptr++) = nodes.size();

   typename NodeContainer::const_iterator i_nodes;
   for ( i_nodes=nodes.begin(); i_nodes!=nodes.end(); ++i_nodes ) {
      (*i_nodes).putToIntBuffer(ptr);
      ptr += node_com_buf_size;
   }

   return;
}


/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
const int *LayerNodeSet<DIM>::acquireNonlocalNodes_unpack( const int *recv_mesg,
                                                            const int src_rank )
{
   /*
    * Unpack node info from recv_mesg into d_nodes,
    * starting at the offset location.
    */
   int node_com_buf_size = LayerNode<DIM>::commBufferSize();

   NodeContainer &nodes = d_nodes[src_rank];
   if ( nodes.size() != 0 ) {
      TBOX_ERROR("Attempting to add nonlocal node when they\n"
                 << "already exist.  Library error!\n");
   }

   const int *ptr = recv_mesg;
   const int n_self_nodes = *(ptr++);

   int i;
   LayerNode<DIM> node;

   for ( i=0; i<n_self_nodes; ++i ) {
      node.getFromIntBuffer(ptr);
      nodes.insert(nodes.end(), node);
      ptr += node_com_buf_size;
   }

   return ptr;
}


/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::acquireNonlocalNodes_unpack(
   const tbox::Array<int> &recv_mesg,
   tbox::Array<int> &proc_offset )
{
   /*
    * Unpack node info from recv_mesg into d_nodes,
    * starting at the offset location.
    * Advance the proc_offset past the used data.
    */
   const int nproc = tbox::SAMRAI_MPI::getNodes();
   int n;
   int node_com_buf_size = LayerNode<DIM>::commBufferSize();

   for ( n=0; n<nproc; ++n ) {
      if ( n != d_rank ) {

         NodeContainer &nodes = d_nodes[n];
         if ( nodes.size() != 0 ) {
            TBOX_ERROR("Attempting to add nonlocal node when they\n"
                       << "already exist.  Library error!\n");
         }

         const int *ptr = recv_mesg.getPointer() + proc_offset[n];
         const int n_self_nodes = *(ptr++);
         proc_offset[d_rank] += (n_self_nodes) * node_com_buf_size;

         int i;
         LayerNode<DIM> node;

         for ( i=0; i<n_self_nodes; ++i ) {
            node.getFromIntBuffer(ptr);
            nodes.insert(nodes.end(), node);
            ptr += node_com_buf_size;
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
void LayerNodeSet<DIM>::setTo( const hier::PatchLevel<DIM> &level )
{
   getMutableNodeContainer().clear();
   if ( d_parallel_state != DISTRIBUTED ) {
      setParallelState(DISTRIBUTED);
   }
   d_ratio = level.getRatio();
   NodeContainer &node_container = getMutableNodeContainer();
   for ( typename hier::PatchLevel<DIM>::Iterator p(level); p; p++ ) {
      const hier::Patch<DIM> &patch = *level.getPatch(*p);
      const hier::Box<DIM> &box = patch.getBox();
      LayerNode<DIM> node( box, *p, d_rank );
      node_container.insert( node_container.end(), node );
   }
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
typename LayerNodeSet<DIM>::NodeContainer::iterator
LayerNodeSet<DIM>::addBox( const hier::Box<DIM> &box,
                           const bool use_vacant_index )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_parallel_state != DISTRIBUTED ) {
      TBOX_ERROR("Individually adding nodes is a local process\n"
                 <<"so it can only be performed in\n"
                 <<"distributed state.");
   }
#endif

   NodeContainer &node_container = getMutableNodeContainer();

   if ( node_container.size() == 0 ) {
      Node new_node = Node( box, 0, tbox::SAMRAI_MPI::getRank() );
      return (node_container.insert(new_node)).first;
   }
   else {
      typename NodeContainer::iterator ni = node_container.end();
      --ni;
      typename Node::LocalIndex new_index = (*ni).getLocalIndex() + 1;
      if ( use_vacant_index ) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( new_index >= 0 );
#endif
         if ( size_t(new_index) != node_container.size() ) {
            /*
             * There is a smaller unused index we can use for the new index.
             */
            for ( new_index=0, ni=node_container.begin();
                  ni!=node_container.end();
                  new_index++, ++ni ) {
               if ( new_index != (*ni).getLocalIndex() ) {
                  break;
               }
            }
            // We should have found an unused index.
            TBOX_ASSERT( ni != node_container.end() );
         }
      }
      Node new_node = Node( box, new_index, d_rank );
      return node_container.insert(ni, new_node);
   }
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void
LayerNodeSet<DIM>::eraseNode( const typename NodeContainer::iterator &inode )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_parallel_state != DISTRIBUTED ) {
      TBOX_ERROR("Individually erasing nodes is a local process\n"
                 <<"so it can only be performed in\n"
                 <<"distributed state.");
   }
#endif

   NodeContainer &node_container = getMutableNodeContainer();
   node_container.erase(inode);
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void
LayerNodeSet<DIM>::eraseNode( const Node &node )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_parallel_state != DISTRIBUTED ) {
      TBOX_ERROR("Individually erasing nodes is a local process\n"
                 <<"so it can only be performed in\n"
                 <<"distributed state.");
   }
#endif

   NodeContainer &node_container = getMutableNodeContainer();
   typename NodeContainer::iterator inode = node_container.find(node);
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( inode == node_container.end() ) {
      TBOX_ERROR("Node to be erased is NOT a part of the node layer.\n");
   }
#endif
   node_container.erase(inode);
   return;
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
bool LayerNodeSet<DIM>::hasNode( const int local_index,
                                 const int owner ) const
{
   hier::Box<DIM> box;
   Node node(box, local_index, owner<0?d_rank:owner);
   return hasNode(node);
}

/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
bool LayerNodeSet<DIM>::hasNode( Node &node ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_parallel_state == DISTRIBUTED && node.getOwnerRank() != d_rank ) {
      TBOX_ERROR("LayerNodeSet: Cannot find remote node while in\n"
                 <<"distributed mode.  LayerNodeSet must be in\n"
                 <<"collected mode first.");
   }
#endif
   const NodeContainer &node_container = getNodeContainer( node.getOwnerRank() );
   typename NodeContainer::const_iterator ni = node_container.find(node);
   if ( ni == node_container.end() ) {
      return false;
   }
   else {
      node.getBox() = (*ni).getBox();
      return true;
   }
}





/*
***********************************************************************
***********************************************************************
*/
template<int DIM>
void LayerNodeSet<DIM>::setTo( const LayerNodeSet &layer )
{
   d_ratio = layer.d_ratio;

   if ( d_parallel_state == DISTRIBUTED ) {
      /*
       * Set only local nodes.
       * Clear neighbor data, which is made obsolete by
       * resetting the nodes.
       */
      NodeContainer &nodes = getMutableNodeContainer();
      nodes.clear();
      nodes = layer.getNodeContainer();
   }

   else /* d_parallel_state == GLOBALIZED */ {
      if ( layer.getParallelState() == DISTRIBUTED ) {
         /*
           Temporarily set state to DISTRIBUTED for the copy.
           Then change the state back.
         */
         deallocateData();
         setParallelState(DISTRIBUTED);
         NodeContainer &nodes = getMutableNodeContainer();
         nodes = layer.getNodeContainer();
         setParallelState(GLOBALIZED);
      }
      else {
         const int nproc = tbox::SAMRAI_MPI::getNodes();
         int n;
         for ( n=0; n<nproc; ++n ) {
            NodeContainer &nodes = d_nodes[n];
            nodes.clear();
            nodes = layer.d_nodes[n];
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
void LayerNodeSet<DIM>::printClassData( std::ostream &co, int detail_depth ) const
{
   int nprocs = tbox::SAMRAI_MPI::getNodes();
   int n, global_size;
   if ( getParallelState() == GLOBALIZED ) {
      for ( n=0, global_size=0; n<nprocs; ++n )
         global_size += getNodeContainer(n).size();
   }
   else {
      global_size = -1;
   }
   const NodeContainer &nodes = getNodeContainer();
   co << "Parallel state: " << (getParallelState()==DISTRIBUTED?"DIST":"GLOB") << '\n'
      << "Ratio         : " << getRefinementRatio() << '\n'
      << "Node count    : " << nodes.size() << ", " << global_size << '\n'
      ;
   if ( detail_depth > 0 ) {
      typename NodeContainer::const_iterator i_node;
      co << "Nodes:\n";
      if ( getParallelState() == GLOBALIZED ) {
         /*
          * Print nodes from all ranks.
          */
         for ( n=0; n<d_rank; ++n ) {
            const NodeContainer &nl_nodes = getNodeContainer(n);
            for ( i_node=nl_nodes.begin(); i_node!=nl_nodes.end(); ++i_node ) {
               Node node = *i_node;
               co << "\t "
                  << node << "\t"
                  << node.getBox().numberCells() << '\n';
            }
         }
      }
      else {
         /*
          * Print local nodes only.
          */
         for ( i_node=nodes.begin(); i_node!=nodes.end(); ++i_node ) {
            Node node = *i_node;
            co << "    "
               << node << "   "
               << node.getBox().numberCells() << '\n';
         }
      }
   }
   return;
}


}
}
#endif
