/*
  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/dlbg/LayerNodeSet.h $
  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
  Revision:    $LastChangedRevision: 2132 $
  Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
  Description: Set of layer nodes in a distributed box graph.
*/

#ifndef included_hier_LayerNodeSet
#define included_hier_LayerNodeSet

#include "SAMRAI_config.h"

#include "LayerNode.h"

#include "PatchLevel.h"

#include "tbox/Array.h"

#include "tbox/PIO.h"

#include "tbox/Pointer.h"

#include "tbox/Timer.h"

#include "tbox/TimerManager.h"

#include <iostream>
#include <map>
#include <set>

namespace SAMRAI {
namespace hier {

/*!
  @brief Encapsulates a set of LayerNode objects on the same index space.

  This class is a part of the distributed box-graph management.
  The distributed box-graph is described in the LayerEdgeSet
  documentation.

  The primary purpose of this class is to organize local node data
  and (when applicable) nonlocal node data.  It performs the
  communication necessary to acquire data on nonlocal nodes.

  A @b local LayerNode is owned by the local process
  (see LayerNode for ownership).
  A @b remote LayerNode is owned by a remote process.

  An LayerNodeSet can be in one of two parallel states:
  - @b DISTRIBUTED: The local process knows only the local nodes
       in the set.
  - @b GLOBALIZED: All processes know all nodes in the set.

  The parallel state is changed by calling setParallelState().
  Note that going from DISTRIBUTED to GLOBALIZED state requires
  an all-to-all gather communication, the performance of which
  should be carefully considered if used frequently.
  The GLOBALIZED state also requires more memory.
  Going from GLOBALIZED state to distributed state is cheap.

  The general attributes of a LayerNodeSet are
  - the set of unique LayerNode objects,
  - the refinement ratio defining their index space, and
  - the parallel state.
  LayerNode objects uniqueness are based on their equality operator,
  which compares both owners and indices.  Therefore,
  a valid set does not contain two nodes with the same owner
  and index.
*/
template<int DIM> class LayerNodeSet
{

public:

   typedef LayerNode<DIM> Node;
   /*!
     @brief Container for nodes.

     This is a sorted container so it can be compared without
     expensive searches.  A node can be removed or added without
     changing the indices of existing nodes.
   */
   typedef std::set<Node> NodeContainer;

   /*!
     @brief Default constructor.

     The default constructor creates object in distributed state.
   */
   LayerNodeSet();

   /*!
     @brief Construct using a deep copy.

     New object has the same parallel state as original.
   */
   LayerNodeSet(const LayerNodeSet &r);

   /*!
     @brief Destructor.

     Deallocate internal data.
   */
   virtual ~LayerNodeSet(void);



   /*!
     @brief Names of parallel states.
   */
   enum ParallelState { DISTRIBUTED, GLOBALIZED };



   /*!
     @brief Assignment operator duplicates edge data and reference
     to layer node set and sets up a similar partner relationship
     and parallel state.

     If @c r is attached to itself, self-attach.  If @c r is attached
     but not to itself, create a duplicate partner and attach to it.

     All other data is directly copied.
   */
   LayerNodeSet &operator=( const LayerNodeSet &r );



   //@{
   /*!
     @name Manually access internal data
   */

   /*!
     @brief Returns the container of nodes for a given process.

     If rank is omitted, it refers to the local process.
     If rank is not the local process,
     the object must be if in GLOBALIZED state.
   */
   const NodeContainer &getNodeContainer( const int rank=-1 ) const;



   /*
     @brief Get const access to layer refinement ratio
     (with respect to a reference level).
   */
   const hier::IntVector<DIM> &getRefinementRatio() const;

   /*
     @brief Set refinement ratio defining the index space
     used by the nodes.

     All nodes owned by the object are assumed to correspond
     to the ratio.
   */
   void setRefinementRatio( const hier::IntVector<DIM> &ratio );


   /*
     @brief Set the internal data to correspond to a
     hier::PatchLevel<DIM> object.

     The internal node data and ratios are set to
     correspond to the given level, replacing current data.
     The parallel state is distrbuted because the level is
     not expected to carry all remote data.
   */
   void setTo( const hier::PatchLevel<DIM> &level );

   /*
     @brief Set the internal data to correspond to a
     hier::PatchLevel<DIM> object.

     The internal node data, ratio are copied from the
     given layer, replacing current data.

     If parallel mode remains unchanged.
     If the parallel mode is GLOBALIZED and the input's is DISTRIBUTED,
     it is temporarily set to DISTRIBUTED for the copy and followed
     by a call to change the mode back to GLOBALIZED.
   */
   void setTo( const LayerNodeSet &layer );

   /*!
     @brief Create new node from given box and append.

     The new node will be assigned an unused local index.
     To be efficient, no communication will be used.
     Therefore, the state must be distributed.

     @return iterator to the new node
   */
   typename NodeContainer::iterator addBox( const hier::Box<DIM> &box,
                                            const bool use_vacant_index = true );

   /*!
     @brief Erase an existing node.

     The given iterator @em MUST be a valid iterator pointing to
     a node currently in this object.
   */
   void eraseNode( const typename NodeContainer::iterator &inode );

   /*!
     @brief Erase an existing node.

     The given node @em MUST match a node currently in this object.
   */
   void eraseNode( const Node &node );

   //@}



   /*!
     @brief Whether object has a given node.

     If @c owner is omitted, the local process is the presumed owner.
     If @em node is not locally owned, the state must be GLOBALIZED.

     @return Whether the specified node exists in the object.
   */
   bool hasNode( const int local_index, const int owner = -1 ) const;

   /*!
     @brief Whether object has a given node.

     Search for a node with the owner and local index
     matching the given node.  (The box of the given
     node is ignored and need not be set.)  If @em node
     is not locally owned, and the state
     not GLOBALIZED, it is an error.

     If the node is found, the box of @c node will be
     set to the correct box for the corresponding node
     and the return value is true.  Otherwise, node is
     not modified and false is returned.

     @return Whether the specified node exists in the object.
   */
   bool hasNode( Node &node ) const;



   /*!
     @brief Set the parallel state.

     This method is not necessarily trivial.
     Acquiring nonlocal node information (when going
     to GLOBALIZED mode) triggers all-gather communication.
     More memory is required to store additional nodes.

     Data not used by the new state gets deallocated.
   */
   void setParallelState( const ParallelState parallel_state );

   /*!
     @brief Set the parallel state for multiple layers.

     This method will combine messages for all the layers
     into a single communication phase.

     This method is stateless.
   */
   void setParallelState( LayerNodeSet *layers,
                          const int num_layers,
                          ParallelState parallel_state) const;

   ParallelState getParallelState() const;


   /*!
     @brief Deallocate internal node data.

     The containers are not deallocated.  They are just emptied.
   */
   void deallocateData();


   virtual void printClassData( std::ostream &os, int detail_depth=0 ) const;



private :


   /*!
     @brief Returns the container of nodes for a given process.

     If rank is omitted, it refers to the local process.
     If rank is not the local process,
     the object must be in GLOBALIZED state.

     This method is like the public getMutableNodeContainer,
     but returns a non-const reference for internal use.
   */
   NodeContainer &getMutableNodeContainer( const int rank=-1 );


   //@{

   /*!
     @brief Get info on nonlocal nodes.

     This requires global communication (all gather).
     This method is more complex than just performing a global data-share.
     To combine communication when a partner is attached,
     both communications are done together.
   */
   void acquireNonlocalNodes( const int hogging_process = -1 );

   //! @brief Pack local nodes into an integer array.
   void acquireNonlocalNodes_pack( tbox::Array<int> &send_mesg,
                                   int offset ) const ;

   //! @brief Unpack nodes from an integer array into internal storage.
   void acquireNonlocalNodes_unpack( const tbox::Array<int> &recv_mesg,
                                     tbox::Array<int> &proc_offset );
   void acquireNonlocalNodes_unpack( const tbox::Array<int> &recv_mesg );
   const int *acquireNonlocalNodes_unpack( const int *recv_mesg,
                                           const int src_rank );

   /*!
     @brief Get info on nonlocal nodes for multiple
     LayerNodeSet objects.

     This method combines communication for the multiple layers
     to increase message passing efficiency.

     Note: This method is stateless (could be static).
   */
   void acquireNonlocalNodes( const int num_sets,
                              LayerNodeSet *layer_node_sets[],
                              const int hogging_process = -1 );
   //@}



   /*!
     @brief Locally-stored nodes.

     In distributed mode, this points to the container
     for the local nodes.  In other modes, this is an
     array of containers, one for each process.
   */
   NodeContainer *d_nodes;


   /*!
     @brief Refinement ratio from some reference such as level 0.
   */
   hier::IntVector<DIM> d_ratio;


   /*!
     @brief State flag.

     Modified by setParallelState().
   */
   ParallelState d_parallel_state;


   /*!
     @brief Process rank (id).

     We save this because we use it frequently.
   */
   const int d_rank;

   /*!
     @brief Timer
   */

   tbox::Pointer<tbox::Timer> t_acquire_nonlocal_nodes;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "LayerNodeSet.I"
#endif

#endif  // included_hier_LayerNodeSet

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LayerNodeSet.C"
#endif
