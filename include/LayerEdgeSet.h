/*
  File:        $RCSfile$
  Copyright:   (c) 1997-2005 The Regents of the University of California
  Revision:    $Revision: 346 $
  Modified:    $Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
  Description: Set of edges in distributed box graph.
*/

#ifndef included_hier_LayerEdgeSet
#define included_hier_LayerEdgeSet

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_LayerNode
#include "LayerNode.h"
#endif

#ifndef included_hier_LayerNodeSet
#include "LayerNodeSet.h"
#endif

#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif

#ifndef included_tbox_MPI
#include "tbox/MPI.h"
#endif

#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#include <iostream>
#include <vector>
#include <map>

using namespace std;

namespace SAMRAI {
namespace hier {

#undef REQUIRE_PARTNER_PARALLEL_STATE_MATCH

/*!
  @brief Encapsulates a set of DLBG edges that connect
  two LayerNodeSet objects.

  @b Distributed @b Layered @b Box-graph (DLBG):

  A box-graph is a graph where the @b nodes correspond to
  boxes and @b edges connect neighboring nodes.  Two nodes
  are @b neighbors if their boxes @b overlap when one is
  grown by some nominal @b ghost-cell-width.
  Each edge is defined by the two nodes it connects.
  The box-graph is useful in SAMR because it shows which
  boxes directly interact with each other.

  We can organize nodes of a box-graph into groups of
  boxes defined on the same index space.  When sorted
  from coarsest to finest index spaces, the nodes fall
  into distinct @b layers.
  A "layer" is similar to a level in the hierarchy in that
  patches in a level share the same index space.
  We use layer instead of level to differentiate the graph from
  the hierarchy and because a layer does not necessarily have
  a one-to-one association with a level.

  For convenience, we define
  - @b layer-node: A set of nodes in same layer (defined
    on the same index space).  This is encapsulated in the
    class LayerNodeSet.
  - @b layer-edge: Set of edges incident from a layer-node.
    This is encapsulated in the class LayerEdgeSet.

  An edge is directional, like an arrow, starting at a @b base
  node and ending at a @b head node.  Each pair of neighboring
  nodes has two edges, pointing in opposite directions.
  This goes also for LayerNodeSet and LayerEdgeSet objects.
  A single LayerEdgeSet object represents a set of directed edges
  incident from one base LayerNodeSet to one head LayerNodeSet.
  The edges pointing in the oposite direction are stored in a
  separate LayerEdgeSet object.

  An object of this class can be @e partnered with another object
  that represents the directed edges in the opposite direction.
  The head of one partner is the base of the other.
  The partner does not provide additional information
  but it does organize the information differently.
  If the edge data is represented as a matrix,
  then one partner is the transpose of the other.
  A partnership allows both objects to be modified
  simultaneously so that their data is consistent.
  This avoids having to modify one after the other,
  which is problematic because the second modification
  begins with inconsistent data.

  Two layers may contain edges in opposite directions and
  not be attached; but they would not be automatically
  coordinated in that case.

  In parallel, a LayerNodeSet can also be in distributed
  or globalized (see LayerNodeSet).  A locally-owned node
  is called a local node.  Nodes owned by other processors
  are remote notes.  Remote nodes are stilled called remote
  nodes when the LayerNodeSet is globalized.
  Local knowledge of a remote node does not make it local.

  In parallel, we may partition the layer-edge and store it
  across many processors.  The edge are classified as follows:
  - @b Local @b edges are incident from a local base
    and incident to two local head.
  - @b Semi-local @b edges are incident from a local base
    and incident to a remote head.
  - @b Remote @b edges connect two remote nodes.

  The @b parallel @b state of a LayerEdgeSet object can be:
  -# @b DISTRIBUTED: An object in this state stores only
     edges incident from local nodes.
  -# @b GLOBALIZED: In this state, the edges are globally
        duplicated on all processors.

  The distribution state is changed by setParallelState().
  Going to a more distributed state primarily means deallocating
  data, but going to a more serial state requires sharing
  of data (communicating).  Currently, collected state for
  LayerEdgeSet objects seems to be unneeded.

  The general attributes of a LayerEdgeSet are:
  - the base LayerNodeSet (a reference to an external object).
  - the ghost cell width to grow the base nodes when checking
    overlaps.
  - the connectivity (table of neighbors for nodes in the
    base LayerNodeSet).
  - the attached partner, if any.
  - the parallel state.

  The connectivity depends on the base LayerNodeSet and the
  ghost-cell-width, so it is discarded when either of these
  is manually changed.
*/
template<int DIM>
class LayerEdgeSet
{

public:

   typedef LayerNode<DIM> Node;
   typedef typename Node::LocalIndex LocalIndex;

   /*!
     @brief Container for nodes.

     This is a sorted container so that it can be compared without
     expensive searches.  A node can be removed or added without
     changing the indices of existing nodes.
   */
   typedef typename LayerNodeSet<DIM>::NodeContainer NodeContainer;

   /*!
     @brief Container for neighbors.

     @internal Neighbors are nodes, so the neighbor container is
     the same as the node container.  For flexibility and readability,
     we use different names.
   */
   typedef typename LayerNodeSet<DIM>::NodeContainer NabrContainer;

   /*!
     @brief Connectivity data between two node layers.

     The connectivity @c m maps a node with index @c i to a set
     of nodes @c m[i], usually its neighbors.
     (The neighbors data is a container of nodes.)

     Because Connectivity is indexed by the LocalIndex instead of
     a combined owner and LocalIndex, each Connectivity object
     is implicitly associated with just one owner process.
   */
   typedef map<LocalIndex,NabrContainer> Connectivity;

   /*!
     @brief Constructor.

     The default constructor creates empty object in
     distributed state.
   */
   LayerEdgeSet();


   /*!
     @brief Destructor.

     Deallocate internal data.
   */
   virtual ~LayerEdgeSet(void);



   /*!
     @brief Names of parallel states.
   */
   enum ParallelState { DISTRIBUTED, GLOBALIZED };
   // typedef typename LayerNodeSet<DIM>::ParallelState ParallelState;



   /*
     @brief Returns the container of nodes for a given process.

     If rank is omitted, it refers to the local process.
     If rank is not the local process,
     the layer node set that was given in initialize()
     must be in GLOBALIZED mode.

     @return container of nodes owned by the given process.
   */
   const NodeContainer &getNodeContainer(const int rank=-1) const;

   /*
     @brief Access connectivity data for nodes owned by a given process.

     This is const because it is public.

     @return container of neighbor sets owned by the given process.
   */
   const Connectivity &getConnectivity( const int rank=-1 ) const;

   /*
     @brief Get const access to layer refinement ratio (with
     respect to a reference level).
   */
   const hier::IntVector<DIM> &getBaseRefinementRatio() const;

   /*
     @brief Get const access to layer refinement ratio (with
     respect to a reference level).
   */
   const hier::IntVector<DIM> &getHeadRefinementRatio() const;


   /*!
     @brief Mutually attach to another layer-edge object
     to create a bi-directional set of edges.

     @b Preconditions on the partner and this object:
     -# The base refinement ratio of one partner and the head
        ratio of the other must be equal, and vice versa.
     -# The product of the base refinement ratio and the
        ghost-cell-width musb be the same for the partners.
        This ensures that the overlap check gives
        the same result regardless of the layer
        whose boxes are grown for the overlap check.
        For the extent of the partnership, this requirement
        must be met.  Resetting the ghost-cell-width of one
        partner forces the ghost-cell-width of the partner
        to change to meet this requirement.

     Both objects must be in the same parallel state.
     Their parallel states will change together for the
     the duration of the partnership.

     Both objects are modified to reference each other as
     mutual partners.
     Both objects must be unpartnered.
     Multiple partners are not allowed.
     Once two layer edges objects are
     attached, they are synchronized when either calls
     a synchronizing method.  Synchronizing means making
     all edges between them seen by both objects.  (This
     does not mean that remote edges are seen in DISTRIBUTED
     modes though.)

     The relationship between two attached objects is symmetric,
     that is, @c a.attach(b) has the same effect as @c b.attach(a).
     It is permissible to pass *this in for partner, which is how
     peer edges are set up.

     Since the partner is not created internally, it is NOT
     explicitly deleted when this object goes out of scope.
     It will simply be detached.
     This is in contrast to createPartner().
   */
   void attachPartner( LayerEdgeSet &partner );


   /*
     @brief Create a partner and attach to it.

     The created partner uses @c partner_base
     for its base layer node set, if given, as though
     initialize() was called with the same argument.

     The created partner satisfies all the criteria
     stated for attachPartner().

     Since the partner is created internally,
     it is deleted when the object goes out of scope.
     This is in contrast to attachPartner().

     See attachPartner()

     If @c connectivity is given, it is assumed to be
     the correct connectivity for the partner.
     If not given, you can use findEdges() to compute
     it after using this method.

     @return the created partner.
   */
   LayerEdgeSet &createPartner(
      const LayerNodeSet<DIM> &partner_base,
      const Connectivity *partner_connectivity = NULL );


   /*!
     @brief Access the attached partner.

     Error if no partner is attached.
   */
   LayerEdgeSet &getPartner();

   /*!
     @brief Access the attached partner.

     Error if no partner is attached.
   */
   const LayerEdgeSet &getPartner() const;

   /*!
     @brief Set data defining the edge set.

     The object must not be attached to a partner.
     This method resets enough data to make coordination
     with partners irrelevant.
     You must detach the partner (if any) before calling this method.
     There should be no reason to use this method while attached.

     A reference is made to the given LayerNodeSet.
     The ghost cell width is copied.
     If @c connectivity is given, a copy is made from it.
     This is meant for use when you have an external means
     for computing the connectivity data.  (The accuracy
     and consistency of the connectivity is not checked!)
     If @c connectivity is not given, it should be created
     using findEdges() before calling a method that uses it.

     The parallel state is initialized to the given @c parallel_state.
     If state is GLOBALIZED:
      - @c base_layer must be in GLOBALIZED mode.
      - connectivity (if given) must be an array the length
        of the number of processes, describing the
        local and semi-local edges for each process.

     A reference is made to the given LayerNodeSet,
     replacing the current node reference.

     When debug is enabled,
     some checks are made to assert that base_layer and
     connectivity are consistent with each other.
   */
   void initialize(
      ParallelState parallel_state,
      const LayerNodeSet<DIM> &base_layer,
      const hier::IntVector<DIM> &head_refinement_ratio,
      const hier::IntVector<DIM> &gcw = hier::IntVector<DIM>(1),
      const Connectivity *connectivity = NULL );


   /*!
     @brief Assignment operator duplicates edge data and reference
     to layer node set and sets up a similar partner relationship
     and parallel mode.

     If @c r is attached to itself, self-attach.  If @c r is attached
     but not to itself, create a duplicate partner and attach to it.

     All other data is directly copied.
   */
   LayerEdgeSet &operator=( const LayerEdgeSet &r );


   //@{
   /*!
     @name Algorithms for computing edge data
   */

   /*!
     @brief Discover and add edges between self and partner.

     Obviously, a partner must be attached.
     Neighbor information on both the object and its partner is modified.
     The effect is the same as if findEdges() were called for the partner.

     This method simply calls findEdges(LayerEdgeSet &a, LayerEdgeSet &b)
     with the partner in the argument.  See that method for communication
     implications.

     If the partner is the same object (finding peer edges),
     edges between a node and itself will be ignored.
   */
   void findEdges();

   /*
     @brief Discover and add edges between two layers.

     The set of edges found depends on the parallel state
     (see setParallelState()).  If in GLOBALIZED mode,
     the entire graph will be found.  If in hybrid mode,
     only the local and semilocal parts will be found.

     The input layers need not be partners, but they
     should be initialized so that they can be attached,
     i.e., the product of the base refinement ratio
     and ghost cell width should be the same for both.

     GLOBALIZED base layer-nodes of the inputs are
     used in finding the edges.

     Edges found are added to appropriate neighbor lists.
     No edge is removed.  If existing edges are invalid,
     remove them first.

     If the node layers of the object and @c other are not in
     hybrid mode, temporay node layers in hybrid mode are created,
     requiring communication.  It is more efficient if they
     are already in hybrid mode.

     If the bases of a and b are the same (finding peer edges),
     edges between a node and itself will be ignored.
   */
   static void findEdges( LayerEdgeSet &a, LayerEdgeSet &b );

   /*!
     @brief Create the layer edge data by bridging
     the heads of two layer edges whose base nodes are identical
     to each other.

     @verbatim

        (bridge base)  ------------bridge------------->  (bridge head)
                       <-------reverse bridge----------        
              ^                                                ^
              |                                                |
              |                                                |
              |                                                |
              +--edge to base-- (middle layer) --edge to head--+

     @endverbatim

     The bases of the two given layers are known as the middle layer.
     For each node in this layer, the two given layers give the sets
     of neighboring nodes at the head and base layers of the current
     object.  These two sets are checked against each other for
     node intersections.

     Conditions:
     -# @c edge_to_head and @c edge_to_base must be incident from
        identical sets of nodes.
     -# @c edge_to_head and @c edge_to_base must have attached partners.
     -# Any parallel state is allowed, but states must be consistent
        between all LayerEdgeSet objects.
     -# At least one of the head or base layers must be completely
        nested in the middle layer.  If neither is nested, some edges
        may be missed.

     To find peer edges, use identical layer-edges for the leg
     to head and leg to base.

     A partner must be attached.
   */
   void bridge( const LayerEdgeSet &edge_to_head,
                const LayerEdgeSet &edge_to_base );


   /*!
     @brief Replace some nodes with others.

     Given the LayerEdgeSet object @c map,
     which maps some nodes in the head of the current
     LayerEdgeSet to new nodes (also in the head),
     replace each given node in the head layer by
     its neighbors in the tail of @c map.

     @verbatim
                               map
                 (old head) ---------> (new head)
                          ^            ^
                           \          /
     edge before mapping -> \        / <- edge after mapping
                             \      /
                              \    /
                              (base)
     @endverbatim

     Representing the map by a LayerEdgeSet
     implicitly requires that any set of nodes
     replacing an old node must overlap the old ones they replace.
     A special algorithm taking advantage of this restriction
     is used to update this object.  If you want to manipulate
     the edge data generally, @em DO @em NOT use this method.

     @c map is a layer edge from nodes to be removed to nodes
     to be in its place.  If any edges in @c map references
     nodes not locally owned, communication with the owner
     is triggered.  In these cases, it is @em CRITICAL that
     the same mapping is introduced on the owner processes.

     We could do the same with a bridge() operation
     centered at the tail of the map.
     but that requires that @c map maps the unchanged
     nodes in the tail to itself, and it also requires
     doing work to map unchanged nodes to themselves.
     This is a lot of extra overhead when map is very small.
     If @c map maps each unchanged node to itself,
     @c a.replaceNodes(map) is equivalent to @c a.bridge(map,a.getPartner()).

     The partner's base will be replaced by the new head.

     Conditions:
     -# @c map must be based at the head of the current object.
     -# @c map and the current object must have attached partners.
     -# Any parallel state is allowed, but states must be consistent
        between all LayerEdgeSet objects.
     -# @c map must not map any node outside the range of the node's
           current box.

     @param map Mapping between old head and new head.
                A partner must be attached.
   */
   void replaceNodes( const LayerEdgeSet &map );

   //@}



   /*!
     @brief Set the parallel distribution state.

     Before putting a LayerEdgeSet in a GLOBALIZED state,
     The LayerNodeSet given in initialize()
     must already be in GLOBALIZED mode.
     The layer-node set must always be in GLOBALIZED
     mode as long as the layer-edge set is.

     In GLOBALIZED state,
     the entire layer is independently duplicated on every process.
     (For consistency, it is imperative that the layer
     is operated on uniformly across all processes.)
     All edges (not just local and semilocal) are found
     when operations to find them are executed.

     This method is not necessarily trivial.
     More memory is required to store additional parts of the graph.

     If a partner is attached, the parallel state of the partner is
     affected the same way.

     Changing parallel state preserves existing neighbor information
     used by the new parallel state.  Neighbor information not used
     by the new parallel state gets deallocated.

     For efficiency reasons, this method NEVER implicitly discovers
     new neighbor information after going to a GLOBALIZED state.
     To get this data, call findEdges() explicitly.

     For serial (one processor) runs, there is no difference between
     the parallel states (except for the names), and there is no real
     cost for switching parallel states.
   */
   void setParallelState( const ParallelState parallel_state );

   ParallelState getParallelState() const;


   /*!
     @brief Deallocate internal data on nodes and neighbors.

     If a partner exists, its internal data is deallocated also.
     This method does @em NOT detach the partner.
   */
   void deallocateData();

   /*!
     @brief Set the max ghost cell width used to check for
     box intersections.

     The default max ghost cell width is 1.
     The max ghost cell width affects the overlap comparison,
     so the old edge data is no longer valid.  The old data
     is discarded.  For efficiency reasons, new edge data is
     not automatically generated.  To get this data, call findEdges()
     explicitly.

     For attached LayerEdgeSet objects,
     the ratio of the ghost-cell-width of the two partners
     must be equal to the inverse ratio of their base layer
     refinement ratios.  This ensures that the overlap
     check gives the same result regardless of the layer
     whose boxes are grown for the overlap check.
     If a partner is attached, the partner's ghost-cell-width
     is automatically changed to satisfy this requirement.
   */
   void setMaxGhostCellWidth( const hier::IntVector<DIM> &gcw );

   const hier::IntVector<DIM> &getMaxGhostCellWidth() const;


   //@{
   /*!
     @name For outputs, error checking and debugging.
   */

   void printClassData( ostream &co, int debase_depth=0 ) const;
   void printEdgeStats( ostream &co ) const;

   /*!
     @brief Check that the node and neighbor containers have
     a one-to-one correspondance.
   */
   void checkNodeNabrCorrespondance() const;

   /*!
     @brief Check that the node referenced by connectivity
     match those in the partner's base LayerNodeSet.

     Obviously, the object must have a partner.
     If the partner's base is not GLOBALIZED,
     a temporary copy is made and globalized for checking,
     triggering communication.
   */
   void checkNodeConsistency() const;

   /*!
     @brief Check that nodes referenced by a given connectivity
     object match nodes in the given node containers.

     @param nodes A node container.
     @param cnect Connectivity data for @c nodes.
     @param head_layer LayerNodeSet describing the head
            with respect to the given connectivity.  This object
            must be in GLOBALIZED state.
   */
   static void checkNodeConsistency(
      const NodeContainer &local_nodes,
      const Connectivity &cnect,
      const LayerNodeSet<DIM> &head_layer );

   /*!
     Assert that connectivity data is correct for the base and head
     (partner's base).

     This is an expensive operation and should only be used for debugging.
     It creates a globalized version of the layer-node data
     (if this data is not already globalized),
     uses the data to compute the edges
     and compares the result to the current edges.

     The object must be partnered and be in hybrid mode
     for checking to work.

     Checking is done as follows:
       - Rebuild the edge containers by finding edges in hybrid mode.
       - Check that the rebuilt containers match the existing containers.

     Currently, the rebuilt edges are rebuilt using findEdges().
     Thus, it may be pointless to use this method as a check for that
     method.
   */
   void checkConnectivity( const LayerNodeSet<DIM> &head ) const;

   //@}


   /*!
     @brief Data sorted by a hash.

     Data indexed by the rank of communication partners are
     assembled in this struct, to be dereferenced in a map.
     This avoids expense of referencing multiple maps
     (just in case dereferencing a map is expensive).
   */
   struct CommunicationStruct {
      vector<int> send_mesg;
      tbox::MPI::request send_rqst;
      int send_size;
      bool send_done;
      vector<int> recv_mesg;
      tbox::MPI::request recv_rqst;
      int recv_size;
      bool recv_done;
      CommunicationStruct() : send_size(-1), send_done(false),
                              recv_size(-1), recv_done(false) {}
   };




private:

   /*
     @brief Access connectivity data for nodes owned by a given process.

     This is NOT be public because it is not const.

     @return container of neighbor sets owned by the given process.
   */
   Connectivity &getConnectivity( const int rank=-1 );

   /*!
     @brief Deallocate internal data but do not touch partner's data.
   */
   void deallocateDataSelf();

   /*!
     @brief Allocate a list of neighbors for every node in a given
     node container.

     For each @c node in the @c node_container, add a corresponding
     neighbor container that coresponds to index @c node.getLocalIndex()
     to @c cnect.
   */
   void allocateConnectivity( const NodeContainer &node_container,
                              Connectivity &cnect ) const;

   /*
     @brief Create a copy of a DISTRIBUTED LayerNodeSet and
     change its state to GLOBALIZED.

     The returned object should be deleted to prevent memory leaks.
   */
   static LayerNodeSet<DIM> *copyAndGlobalize( const LayerNodeSet<DIM> &r );

   /*!
     @brief Discover phase of bridging.

     Given a set of edges from a "middle" layer to a "head" layer and
     a set of edges from the same "middle" layer to a "base" layer,
     discover edges between the head and base.

     - GLOBALIZED state: local, semilocal and remote edges
       are discovered locally.
     - DISTRIBUTED state: local, semilocal and remote edges
       are discovered locally.

     If partner is self (discovering peer edges), a node's intersection
     with itself is are disregarded.

     Edges connecting nonlocal nodes may be discovered.  Such information
     is appended to messages to be sent to the owners of those nodes.
     See bridge_shareEdges().
   */
   void bridge_discoverEdges(
      const LayerEdgeSet &edge_to_head,
      const LayerEdgeSet &edge_to_base,
      map<LocalIndex,CommunicationStruct> &mesg_to_owners,
      const bool ignore_self_overlap );

   /*!
     @brief Share phase of bridging.

     Messages are created by bridge_discoverEdges().
   */
   void bridge_shareEdges( map<LocalIndex,CommunicationStruct> &mesg_to_owners );

   /*!
     @brief Discover and add edges between self and another,
     including remote edges.

     Neighbor information on both the object and its partner is modified.
     Edges found are added to appropriate neighbor lists.  No edge is
     removed.  If existing edges are invalid, remove them first.

     This method is symmetric.
   */
   void findEdges_serial( LayerEdgeSet &other );

   /*!
     @brief Discover and add edges between self and another,
     including remote edges.

     Neighbor information on both the object and its partner is modified.
     Edges found are added to appropriate neighbor lists.  No edge is
     removed.  If existing edges are invalid, remove them first.

     This method is not symmetric.
   */
   void findEdges_quadratic( const LayerEdgeSet &other );

   /*!
     @brief Discover and add edges from base and externally
     provided head layer.

     Edges found are added to appropriate neighbor lists.  No edge is
     removed.  If existing edges are invalid, remove them first.
   */
   void findEdges_rbbt( const LayerNodeSet<DIM> &head,
                        const bool ignore_self_overlap = false );


   //@{ @name Private utilities.

   /*!
     @brief Get a list of owners of neighbors.
   */
   void getNabrOwnerSet( const Connectivity &cnect,
                         set<int> &owners ) const;

   //@}


   /*!
     @brief Max ghost cell width.

     When determining intersection of two boxes, the box in the
     node layer with finer refinement ratio is grown by this
     amount. If both layers are at the same ratio, it does not
     matter which is grown.
   */
   hier::IntVector<DIM> d_gcw;

   /*!
     @brief Head layer refinement ratio.

     The neighbors in d_cnect correspond to this ratio.
     The partner (if any) must have a base layer with this ratio.
     (The base layer refinement ratio is in d_nodes.)
   */
   hier::IntVector<DIM> d_head_ratio;

   /*!
     @brief Reference to the base layer of nodes.

     This pointer references an object that may be shared
     among many LayerEdgeSet.  Therefore, it will not
     be modified by this class.
   */
   const LayerNodeSet<DIM> *d_base;

   /*!
     @brief Neighbor data.

     In serial mode, this is an array of Connectivity objects,
     one for each process.  In non-serial mode, this is a
     single Connectivity object for the local process.
   */
   Connectivity *d_cnect;


   /*!
     @brief State flag.

     Modified by setParallelState().
   */
   ParallelState d_parallel_state;

   LayerEdgeSet *d_partner;

   /*!
     @brief Whether partner is internally created and managed.

     Set by attachPartner() and createPartner().
   */
   bool d_partner_is_managed;

   const int d_rank;
   const int d_nprocs;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "LayerEdgeSet.I"
#endif

#endif // included_hier_LayerEdgeSet

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LayerEdgeSet.C"
#endif
