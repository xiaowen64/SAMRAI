/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/dlbg/LayerNode.I $
 * Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
  : hier::Box<DIM>(),
    d_owner_rank(0),
    d_local_index(-1)
{
   return;
}

template<int DIM> inline
LayerNode<DIM>::LayerNode( const hier::Box<DIM> &box,
                                  const LocalIndex index,
                                  const int owner )
  : hier::Box<DIM>(box),
    d_owner_rank(owner),
    d_local_index(index)
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
   return (hier::Box<DIM>&)(*this);
}


template<int DIM> inline
const hier::Box<DIM> &LayerNode<DIM>::getBox() const {
   return (const hier::Box<DIM>&)(*this);
}


template<int DIM> inline
bool LayerNode<DIM>::operator==( const LayerNode &r ) const {
   bool rval = ( d_owner_rank == r.d_owner_rank ) && ( d_local_index == r.d_local_index );
#if DEBUG_CHECK_ASSERTIONS
   if ( rval && (hier::Box<DIM>&)(*this) != (hier::Box<DIM>&)(r) ) {
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
   buffer[2] = this -> lower()(0);
   buffer[3] = this -> upper()(0);
   if ( DIM >= 2 ) {
      buffer[4] = this -> lower()(1);
      buffer[5] = this -> upper()(1);
   }
   if ( DIM >= 3 ) {
      buffer[6] = this -> lower()(2);
      buffer[7] = this -> upper()(2);
   }
   return;
}


template<int DIM> inline
void LayerNode<DIM>::getFromIntBuffer( const int *buffer )
{
   d_owner_rank = buffer[0];
   d_local_index = buffer[1];
   this -> lower()(0) = buffer[2];
   this -> upper()(0) = buffer[3];
   if ( DIM >= 2 ) {
      this -> lower()(1) = buffer[4];
      this -> upper()(1) = buffer[5];
   }
   if ( DIM >= 3 ) {
      this -> lower()(2) = buffer[6];
      this -> upper()(2) = buffer[7];
   }
   return;
}



#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
