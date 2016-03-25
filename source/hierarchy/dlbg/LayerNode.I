/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/dlbg/LayerNode.I $
 * Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1859 $
 * Modified:    $LastChangedDate: 2008-01-15 07:55:41 -0800 (Tue, 15 Jan 2008) $
 * Description: Node in the distribued box graph.
 */

#include "LayerNode.h"


#ifndef included_tbox_Utilities
#include <tbox/Utilities.h>
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
LayerNode<DIM>::LayerNode()
   :d_owner_rank(0),
    d_local_index(-1),
    d_box()
{
   return;
}

template<int DIM> inline
LayerNode<DIM>::LayerNode( const hier::Box<DIM> &box,
                                  const LocalIndex index,
                                  const int owner )
  :d_owner_rank(owner),
   d_local_index(index),
   d_box(box)
{
   return;
}


template<int DIM> inline
int LayerNode<DIM>::getOwnerRank() const {
   return d_owner_rank;
}


template<int DIM> inline
int LayerNode<DIM>::getLocalIndex() const {
   return d_local_index;
}


template<int DIM> inline
hier::Box<DIM> &LayerNode<DIM>::getBox() {
   return d_box;
}


template<int DIM> inline
const hier::Box<DIM> &LayerNode<DIM>::getBox() const {
   return d_box;
}


template<int DIM> inline
bool LayerNode<DIM>::operator==( const LayerNode &r ) const {
   bool rval = ( d_owner_rank == r.d_owner_rank ) && ( d_local_index == r.d_local_index );
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( rval && d_box != r.getBox() ) {
      TBOX_ERROR("Inconsistent box data in LayerNode object");
   }
#endif
   return rval;
}


template<int DIM> inline
bool LayerNode<DIM>::operator!=( const LayerNode &r ) const {
   bool rval = ( d_owner_rank != r.d_owner_rank ) || ( d_local_index != r.d_local_index );
   return rval;
}


template<int DIM> inline
bool LayerNode<DIM>::operator<( const LayerNode &r ) const {
   return ( d_owner_rank < r.d_owner_rank ) ||
      ( ( d_owner_rank == r.d_owner_rank ) && ( d_local_index < r.d_local_index ) );
}


template<int DIM> inline
bool LayerNode<DIM>::operator>( const LayerNode &r ) const {
   return ( d_owner_rank > r.d_owner_rank ) ||
      ( ( d_owner_rank == r.d_owner_rank ) && ( d_local_index > r.d_local_index ) );
}


template<int DIM> inline
bool LayerNode<DIM>::operator<=( const LayerNode &r ) const {
   return ( d_owner_rank < r.d_owner_rank ) ||
      ( ( d_owner_rank == r.d_owner_rank ) && ( d_local_index <= r.d_local_index ) );
}


template<int DIM> inline
bool LayerNode<DIM>::operator>=( const LayerNode &r ) const {
   return ( d_owner_rank > r.d_owner_rank ) ||
      ( ( d_owner_rank == r.d_owner_rank ) && ( d_local_index >= r.d_local_index ) );
}



template<int DIM> inline
int LayerNode<DIM>::commBufferSize()
{
   return 2 + 2*DIM;
}


template<int DIM> inline
void LayerNode<DIM>::putToIntBuffer( int *buffer ) const
{
   buffer[0] = d_owner_rank;
   buffer[1] = d_local_index;
   buffer[2] = d_box.lower()(0);
   buffer[3] = d_box.upper()(0);
   if ( DIM >= 2 ) {
      buffer[4] = d_box.lower()(1);
      buffer[5] = d_box.upper()(1);
   }
   if ( DIM >= 3 ) {
      buffer[6] = d_box.lower()(2);
      buffer[7] = d_box.upper()(2);
   }
   return;
}


template<int DIM> inline
void LayerNode<DIM>::getFromIntBuffer( const int *buffer )
{
   d_owner_rank = buffer[0];
   d_local_index = buffer[1];
   d_box.lower()(0) = buffer[2];
   d_box.upper()(0) = buffer[3];
   if ( DIM >= 2 ) {
      d_box.lower()(1) = buffer[4];
      d_box.upper()(1) = buffer[5];
   }
   if ( DIM >= 3 ) {
      d_box.lower()(2) = buffer[6];
      d_box.upper()(2) = buffer[7];
   }
   return;
}



#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
