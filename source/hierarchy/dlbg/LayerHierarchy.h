/*
  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/hierarchy/dlbg/LayerHierarchy.h $
  Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
  Revision:    $LastChangedRevision: 2195 $
  Modified:    $LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
  Description: Box graph representing hierarchy.
*/

#ifndef included_hier_LayerHierarchy
#define included_hier_LayerHierarchy

#include "SAMRAI_config.h"

#include "LayerEdgeSet.h"

#include "LayerNodeSet.h"

#include "PatchLevel.h"

#include "tbox/Array.h"

#include "tbox/PIO.h"

#include "tbox/Pointer.h"

namespace SAMRAI {
namespace hier {


/*!
  @brief Distributed nested-level box graph for a patch hierarchy.

  This object maintains the DNBG edges within a hierarchy.

  A "layer" is similar to a level in the hierarchy in that
  the nodes in a layer share the same structured index space,
  just as patches in a level share the same index space.
  We use layer instead of level to differentiate the graph from
  the hierarchy and because a layer does not necessarily have
  a one-to-one association with a level.
  A "layer edge" is the collection of individual DNBG edges
  that are incident from nodes in the same layer in the DNBG.

*/
template <int DIM>
class LayerHierarchy
{

public:

   typedef LayerEdgeSet<DIM> LayerEdges;
   typedef LayerNodeSet<DIM> LayerNodes;

   /*!
     @brief Constructor.

     The default constructor creates object with distributed seriality.
   */
   LayerHierarchy();

   /*!
     @brief Destructor.

     Deallocate internal data.
   */
   virtual ~LayerHierarchy(void);



   /*
     @brief Get const access to the node container of a level.

     Direct manipulation of nodes is not allowed.
     There are specific methods for doing that.

     @return the layer of nodes for the given level.
   */
   const LayerNodes &getNodes( int level_number ) const;

   /*
     @brief Get const access to peer edges for a level.

     Direct manipulation of neighbor data is not allowed.
     There are specific methods for doing that.

     @return the layer of peer edges for the given level.
   */
   const LayerEdges &getPeerEdges( int level_number ) const;

   /*
     @brief Get const access to fine edges for a level.

     Direct manipulation of neighbor data is not allowed.
     There are specific methods for doing that.

     @return the layer of fine edges for the given level.
   */
   const LayerEdges &getFineEdges( int level_number ) const;

   /*
     @brief Get const access to coarse edges for a level.

     Direct manipulation of neighbor data is not allowed.
     There are specific methods for doing that.

     @return the layer of coarse edges for the given level.
   */
   const LayerEdges &getCoarseEdges( int level_number ) const;



private:

   int d_num_levels;

   tbox::Array<LayerEdgeSet<DIM> > d_peer_edges;
   tbox::Array<LayerEdgeSet<DIM> > d_fine_edges;
   tbox::Array<LayerEdgeSet<DIM> > d_coarse_edges;

   tbox::Array<LayerNodeSet<DIM> > d_node_sets;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "LayerHierarchy.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LayerHierarchy.C"
#endif

#endif  // included_hier_LayerHierarchy
