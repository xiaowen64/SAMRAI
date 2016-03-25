/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/dlbg/LayerNode.h $
 * Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2132 $
 * Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
 * Description: Node in the distribued box graph.
 */

#ifndef included_hier_LayerNode
#define included_hier_LayerNode

#include <iostream>

#include "SAMRAI_config.h"

#include "Box.h"

namespace SAMRAI {
   namespace hier {


/*!
 * @brief Encapsulates the node on a DLBG.
 *
 * A DLBG node is basically a box with a specific owner process
 * and an index on that process.
 *
 * Comparison operators are implemented for sorting nodes
 * and instantiating (STL) sets of LayerNode's.
 * The owners and local indices are used for all comparisons
 * The owners ranks are compared first, followed by the local indices.
 * Less-than and greater-than comparisons are primarily used for
 * sorting nodes.
 */
template<int DIM> class LayerNode
{

public:

   typedef int LocalIndex;

   /*!
    * @brief Constructor.
    */
   LayerNode();

   /*!
    * @brief Constructor.
    */
   LayerNode( const hier::Box<DIM> &box,
                    const LocalIndex index=-1,
                    const int owner_rank=-1 );

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   virtual ~LayerNode(void);


   int getOwnerRank() const;
   LocalIndex getLocalIndex() const;
   hier::Box<DIM> &getBox();
   const hier::Box<DIM> &getBox() const;


   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    *
    * For the equality operator,
    * The box should be the same if data is consistent.
    * If debug is turned on and the boxes do not match
    * while the owners and indices match, then an unrecoverable
    * exception is thrown.
    */
   bool operator==( const LayerNode &r ) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const LayerNode&);
    */
   bool operator!=( const LayerNode &r ) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const LayerNode&);
    */
   bool operator<( const LayerNode &r ) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const LayerNode&);
    */
   bool operator>( const LayerNode &r ) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const LayerNode&);
    */
   bool operator<=( const LayerNode &r ) const;

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const LayerNode&);
    */
   bool operator>=( const LayerNode &r ) const;


   //@}


   //@{
   //! @name Support for message passing
   /*!
    * @brief Give number of ints required in message passing buffer.
    *
    * This number is independent of instance (but dependent on
    * dimension).
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   static int commBufferSize();
   /*!
    * @brief Put self into a int buffer.
    *
    * Number of ints written is given by communicationSize().
    */
   void putToIntBuffer( int *buffer ) const;
   /*!
    * @brief Set self according to data in int buffer.
    *
    * Number of ints read is given by communicationSize().
    */
   void getFromIntBuffer( const int *buffer );
   //@}


   template <int D>
   friend std::ostream &operator<<( std::ostream &co, const LayerNode<D> &r );


private:

   /*!
    * @brief Rank of owner of this node.
    *
    * The DLBG is inherently distributed, so a node always has an owner.
    * One thing that requires the owner be a state variable is that
    * neighbor containers do not explicitly keep track of the owners
    * of the neighbor nodes.
    */
   int d_owner_rank;

   /*!
    * @brief Local index on the owner process.
    *
    * The node may be referenced on processes that do not own it.
    * For example, a process that owns the node's neighbor needs
    * to reference the node.  When communicating with the owner
    * about the node, other processes need to state the local index
    * of the node so the owner knows what node is being talked about.
    */
   LocalIndex d_local_index;

   
   Box<DIM> d_box;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "LayerNode.I"
#endif

#endif  // included_hier_LayerNode

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LayerNode.C"
#endif
