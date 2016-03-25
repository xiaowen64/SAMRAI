//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outernode/OuternodeData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outernode centered patch data type
//

#ifndef included_pdat_OuternodeData_C
#define included_pdat_OuternodeData_C

#include "OuternodeData.h"

#include "Box.h"
#include "BoxList.h"
#include "NodeData.h"
#include "NodeGeometry.h"
#include "NodeOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include <stdio.h>

#define PDAT_OUTERNODEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "OuternodeData.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for outernode data objects.  The           *
* constructor simply initializes data variables and sets up the         *
* array data.                                                           *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuternodeData<DIM,TYPE>::OuternodeData(
   const hier::Box<DIM>& box,
   int depth,
   tbox::Pointer<tbox::Arena> pool)
:  hier::PatchData<DIM>(box, hier::IntVector<DIM>(0)),
   d_depth(depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   for (int d = 0; d < DIM; d++) {

      hier::Box<DIM> nodebox = NodeGeometry<DIM>::toNodeBox(box);

      for ( int dh = d+1; dh < DIM; dh++ ) {

         /*
	  * For dimensions higher than d, narrow the box down to avoid
	  * representing edge and corner nodes multiple times.
          *
          *  i.e.    Y--Y--Y  outernodeX0 defined on nodes (0,1)
          *          |  |  |  outernodeX1 defined on nodes (2,1)
          *          X--o--X  outernodeY0 defined on node  (0,0)-(2,0)
          *          |  |  |  outernodeY1 defined on node  (0,2)-(2,2)
          *          Y--Y--Y
          *         node box
          */
         nodebox.lower(dh) += 1;
         nodebox.upper(dh) -= 1;
      }

      hier::Box<DIM> outernodebox = nodebox;
      outernodebox.upper(d) = nodebox.lower(d);
      outernodebox.lower(d) = nodebox.lower(d);
      if ( outernodebox.size() > 0 ) {
         d_data[d][0].initializeArray(outernodebox, depth, pool);
      }
      outernodebox = nodebox;
      outernodebox.lower(d) = nodebox.upper(d);
      outernodebox.upper(d) = nodebox.upper(d);
      d_data[d][1].initializeArray(outernodebox, depth, pool);

   }
}

template <int DIM, class TYPE>
OuternodeData<DIM,TYPE>::~OuternodeData()
{
}

/*
*************************************************************************
*                                                                       *
* The following are private and cannot be used, but they are defined    *
* here for compilers that require that every template declaration have  *
* a definition (a stupid requirement, if you ask me).                   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuternodeData<DIM,TYPE>::OuternodeData(
   const OuternodeData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::operator=(const OuternodeData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy between an outernode patch data type (source) and *
* a node patch data type (destination) where the index spaces overlap.  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const NodeData<DIM,TYPE>* const t_node_src =
      dynamic_cast<const NodeData<DIM,TYPE> *>(&src);
   const OuternodeData<DIM,TYPE>* const t_onode_src =
      dynamic_cast<const OuternodeData<DIM,TYPE> *>(&src);

   if ( t_node_src != NULL ) {
      copyFromNode( *t_node_src );
   } else if ( t_onode_src != NULL ) {
      copyFromOuternode( *t_onode_src );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy error!\n"
                 << "Can copy only from NodeData<DIM,TYPE> or "
                 << "OuternodeData<DIM,TYPE> of the same "
		 << "DIM and TYPE.");
   }

}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   NodeData<DIM,TYPE> *t_node_dst =
      dynamic_cast<NodeData<DIM,TYPE> *>(&dst);
   OuternodeData<DIM,TYPE> *t_onode_dst =
      dynamic_cast<OuternodeData<DIM,TYPE> *>(&dst);

   if ( t_node_dst != NULL ) {
      copyToNode( *t_node_dst );
   } else if ( t_onode_dst != NULL ) {
      copyToOuternode( *t_onode_dst );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy2 error!\n"
                 << "Can copy only to NodeData<DIM,TYPE> or "
                 << "OuternodeData<DIM,TYPE> of the same "
		 << "DIM and TYPE.");
   }
}

/*
*************************************************************************
*                                                                       *
* Copy data from the source into the destination according to the       *
* overlap descriptor.                                                   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                   const hier::BoxOverlap<DIM>& overlap)
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const NodeData<DIM,TYPE> *t_node_src = 
      dynamic_cast<const NodeData<DIM,TYPE> *>(&src);
   const OuternodeData<DIM,TYPE> *t_onode_src = 
      dynamic_cast<const OuternodeData<DIM,TYPE> *>(&src);

   if ( t_node_src != NULL ) {
      copyFromNode( *t_node_src, *t_overlap );
   } else if ( t_onode_src != NULL ) {
      copyFromOuternode( *t_onode_src, *t_overlap );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy error!\n"
                 << "Can copy only from NodeData<DIM,TYPE> or "
                 << "OuternodeData<DIM,TYPE> of the same "
		 << "DIM and TYPE.");
   }

}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   NodeData<DIM,TYPE> *t_node_dst = 
      dynamic_cast<NodeData<DIM,TYPE> *>(&dst);
   OuternodeData<DIM,TYPE> *t_onode_dst = 
      dynamic_cast<OuternodeData<DIM,TYPE> *>(&dst);

   if ( t_node_dst != NULL ) {
      copyToNode( *t_node_dst, *t_overlap );
   } else if ( t_onode_dst != NULL ) {
      copyToOuternode( *t_onode_dst, *t_overlap );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy2 error!\n"
                 << "Can copy only to NodeData<DIM,TYPE> or "
                 << "OuternodeData<DIM,TYPE> of the same "
		 << "DIM and TYPE.");
   }

}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy from a node data object to this outernode data    *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyDepth(int dst_depth,
                                        const NodeData<DIM,TYPE>& src,
                                        int src_depth)
{
   const ArrayData<DIM,TYPE>& node_array = src.getArrayData();
   for (int d = 0; d < DIM; d++ ) {
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE>& onode_array = d_data[d][loc];
         onode_array.copyDepth(dst_depth,
                               node_array,
                               src_depth,
                               onode_array.getBox());
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy to a node data object from this outernode data    *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyDepth2(int dst_depth,
                                         NodeData<DIM,TYPE>& dst,
                                         int src_depth) const
{
   ArrayData<DIM,TYPE>& node_array = dst.getArrayData();
   for (int d = 0; d < DIM; d++ ) {
      for ( int loc = 0; loc < 2; loc++ ) {
         const ArrayData<DIM,TYPE>& onode_array = d_data[d][loc];
         node_array.copyDepth(dst_depth,
                              onode_array,
                              src_depth,
                              onode_array.getBox());
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Add source data to the destination according to overlap descriptor.   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::sum(
   const hier::PatchData<DIM>& src,
   const hier::BoxOverlap<DIM>& overlap)
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const OuternodeData<DIM,TYPE> *t_onode_src = 
      dynamic_cast<const OuternodeData<DIM,TYPE> *>(&src);

   // NOTE:  We assume this operation is only needed to
   //        copy and add data to another outernode data
   //        object.  If we ever need to provide this for node
   //        data or other flavors of the copy operation, we
   //        should refactor the routine similar to the way 
   //        the regular copy operations are implemented.
   if ( t_onode_src == NULL ) {
      TBOX_ERROR("OuternodeData<DIM>::sum error!\n"
                 << "Can copy and add only from OuternodeData<DIM,TYPE> "
                 << "of the same DIM and TYPE.");
   } else {

      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

      for ( int src_d = 0; src_d < DIM; src_d++ ) {
         for ( int src_p = 0; src_p < 2; src_p++ ) {
 
            const ArrayData<DIM,TYPE> &src_array = 
               t_onode_src->d_data[src_d][src_p];
            const hier::BoxList<DIM>& box_list = 
               t_overlap->getDestinationBoxList();
 
            for ( int dst_d = 0; dst_d < DIM; dst_d++ ) {
               for ( int dst_p = 0; dst_p < 2; dst_p++ ) {
                  if (d_data[dst_d][dst_p].isInitialized()) {
                     d_data[dst_d][dst_p].sum( 
                        src_array, box_list, src_offset);
                  }
               }
            }

         }
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Calculate the buffer space needed to pack/unpack messages on the box  *
* region using the overlap descriptor.                                  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
bool OuternodeData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template <int DIM, class TYPE>
int OuternodeData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   int size = 0;
   const hier::BoxList<DIM>& boxlist = t_overlap->getDestinationBoxList();
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      size += d_data[d][0].getDataStreamSize(boxlist, src_offset);
      size += d_data[d][1].getDataStreamSize(boxlist, src_offset);
   }
   return(size);
}

/*
*************************************************************************
*                                                                       *
* Pack/unpack data into/out of the message streams using the index      *
* space in the overlap descriptor.                                      *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::BoxList<DIM>& dst_boxes    = t_overlap->getDestinationBoxList();
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes); 
           dst_box; dst_box++) {
         const hier::Box<DIM> src_box = hier::Box<DIM>::shift(dst_box(), 
                                                              -src_offset);
         for (int loc = 0; loc < 2; loc++) {
            const hier::Box<DIM> intersect = src_box * d_data[d][loc].getBox();
            if (!intersect.empty()) {
               d_data[d][loc].packStream(stream, 
                                         hier::Box<DIM>::shift(intersect, 
                                                               src_offset), 
                                         src_offset);
            }
         }
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::BoxList<DIM>& dst_boxes = t_overlap->getDestinationBoxList();
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes); 
           dst_box; dst_box++) {
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect = 
               dst_box() * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStream(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Unpack data from the message stream and add to this outernode data    *
* object using the index space in the overlap descriptor.               *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::unpackStreamAndSum(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::BoxList<DIM>& dst_boxes = t_overlap->getDestinationBoxList();
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes);
           dst_box; dst_box++) {
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect =
               dst_box() * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStreamAndSum(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Calculate the amount of memory space needed to represent the data     *
* for a  outernode centered grid.                                       *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
size_t OuternodeData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   size_t size = 0;
   for (int d = 0; d < DIM; d++) {
      hier::Box<DIM> loc0 = NodeGeometry<DIM>::toNodeBox(box);
      hier::Box<DIM> loc1 = NodeGeometry<DIM>::toNodeBox(box);
      loc0.upper(d) = box.lower(d);
      loc0.lower(d) = box.lower(d);
      loc1.lower(d) = box.upper(d);
      loc1.upper(d) = box.upper(d);

      for ( int dh = d+1; dh < DIM; dh++ ) {

         /*
	  * For dimensions higher than d, narrow the box down to avoid
	  * representing edge and corner nodes multiple times.
          */
         loc0.lower(dh) += 1;
         loc0.upper(dh) -= 1;
         loc1.lower(dh) += 1;
         loc1.upper(dh) -= 1;
      }
      size += ArrayData<DIM,TYPE>::getSizeOfData(loc0, depth)
	    + ArrayData<DIM,TYPE>::getSizeOfData(loc1, depth);
   }
   return(size);
}

/*
*************************************************************************
*                                                                       *
* Compute the box of valid node indices given values of                 *
* dimension and side designating the set of data indices.               *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
hier::Box<DIM> 
OuternodeData<DIM,TYPE>::getDataBox(int face_normal, int side)
{
   if ( face_normal < 0 || face_normal >=DIM || side < 0 || side > 1 ) {
      TBOX_ERROR("Bad values for face_normal and/or side in\n"
                 "OuternodeData<DIM>::getDataBox().\n");
   }

   /*
    * We start with the full box and chop it down to the databox
    * corresponding to the given face_normal and side.
    */
   hier::Box<DIM> databox = NodeGeometry<DIM>::toNodeBox(this -> getBox());
   const hier::IntVector<DIM> &ghosts = this -> getGhostCellWidth();

   for ( int dh = face_normal+1; dh < DIM; dh++ ) {

      /*
       * For dimensions higher than d, narrow the box down to avoid
       * representing edge and corner nodes multiple times.
       */
      databox.lower(dh) += 1;
      databox.upper(dh) -= 1;
   }

   if ( side == 0 ) {
      databox.upper(face_normal) = databox.lower(face_normal);
      databox.lower(face_normal) = databox.lower(face_normal) - 
                                   ghosts(face_normal);
   }
   else { // side == 1
      databox.lower(face_normal) = databox.upper(face_normal);
      databox.upper(face_normal) = databox.upper(face_normal) + 
                                   ghosts(face_normal);
   }
   return databox;
}

/*
*************************************************************************
*                                                                       *
* Fill the outernode centered box with the given value.                 *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fill(t, d);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fill(t, d);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::fill(const TYPE& t,
                                   const hier::Box<DIM>& box,
                                   int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fill(t, NodeGeometry<DIM>::toNodeBox(box), d);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fill(t, NodeGeometry<DIM>::toNodeBox(box), d);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::fillAll(const TYPE& t)
{
   for (int i = 0; i < DIM; i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fillAll(t);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fillAll(t);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   for (int i = 0; i < DIM; i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fillAll(t, NodeGeometry<DIM>::toNodeBox(box));
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fillAll(t, NodeGeometry<DIM>::toNodeBox(box));
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy between an outernode patch data type (source) and *
* a node patch data type (destination) where the index spaces overlap.  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyFromNode(const NodeData<DIM,TYPE>& src)
{
   const ArrayData<DIM,TYPE> &node_array = src.getArrayData();
   for ( int d = 0; d < DIM; d++ ) {
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE> &onode_array = d_data[d][loc];
         if (onode_array.isInitialized()) {
            onode_array.copy( node_array, onode_array.getBox() );
         }
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyToNode(NodeData<DIM,TYPE>& dst) const
{
   ArrayData<DIM,TYPE> &node_array = dst.getArrayData();
   for (int d = 0; d < DIM; d++) {
      for ( int loc = 0; loc < 2; loc++ ) {
         if (d_data[d][loc].isInitialized()) {
            node_array.copy(d_data[d][loc], d_data[d][loc].getBox());
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Copy data from the source into the destination according to the       *
* overlap descriptor.                                                   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyFromNode(const NodeData<DIM,TYPE>& src,
                                           const NodeOverlap<DIM>& overlap)
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList();
      if (d_data[d][0].isInitialized()) {
         d_data[d][0].copy( src.getArrayData(), box_list, src_offset);
      }
      if (d_data[d][1].isInitialized()) {
         d_data[d][1].copy( src.getArrayData(), box_list, src_offset);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyToNode(NodeData<DIM, TYPE>& dst,
                                         const NodeOverlap<DIM>& overlap) const
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList();
   for (int d = 0; d < DIM; d++) {
      if (d_data[d][0].isInitialized()) {
         dst.getArrayData().copy(d_data[d][0], box_list, src_offset);
      }
      if (d_data[d][1].isInitialized()) {
         dst.getArrayData().copy(d_data[d][1], box_list, src_offset);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyFromOuternode( 
   const OuternodeData<DIM,TYPE> &src )
{
   for ( int src_d = 0; src_d < DIM; src_d++ ) {
      for ( int src_p = 0; src_p < 2; src_p++ ) {

	 const ArrayData<DIM,TYPE> &src_array = src.d_data[src_d][src_p];

	 for ( int dst_d = 0; dst_d < DIM; dst_d++ ) {
	    for ( int dst_p = 0; dst_p < 2; dst_p++ ) {
               ArrayData<DIM,TYPE> &onode_array = d_data[dst_d][dst_p];
               if (onode_array.isInitialized()) {
                  onode_array.copy( src_array, onode_array.getBox() );
               }
	    }
	 }
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyFromOuternode( 
   const OuternodeData<DIM,TYPE> &src,
   const NodeOverlap<DIM> &overlap )
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   for ( int src_d = 0; src_d < DIM; src_d++ ) {
      for ( int src_p = 0; src_p < 2; src_p++ ) {

	 const ArrayData<DIM,TYPE> &src_array = src.d_data[src_d][src_p];
	 const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList();

	 for ( int dst_d = 0; dst_d < DIM; dst_d++ ) {
	    for ( int dst_p = 0; dst_p < 2; dst_p++ ) {
               if (d_data[dst_d][dst_p].isInitialized()) {
                  d_data[dst_d][dst_p].copy( src_array, box_list, src_offset);
               }
	    }
	 }
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyToOuternode( OuternodeData<DIM,TYPE> &dst ) const
{
   for ( int dst_d = 0; dst_d < DIM; dst_d++ ) {
      for ( int dst_p = 0; dst_p < 2; dst_p++ ) {

	 ArrayData<DIM,TYPE> &dst_array = dst.d_data[dst_d][dst_p];

         for ( int src_d = 0; src_d < DIM; src_d++ ) {
            for ( int src_p = 0; src_p < 2; src_p++ ) {
               if (d_data[src_d][src_p].isInitialized()) {
                  dst_array.copy( d_data[src_d][src_p],
                                  d_data[src_d][src_p].getBox() );
               }
	    }
	 }
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::copyToOuternode( OuternodeData<DIM,TYPE> &dst,
		      const NodeOverlap<DIM> &overlap ) const
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList();

   for ( int dst_d = 0; dst_d < DIM; dst_d++ ) {
      for ( int dst_p = 0; dst_p < 2; dst_p++ ) {

	 ArrayData<DIM,TYPE> &dst_array = dst.d_data[dst_d][dst_p];
         for ( int src_d = 0; src_d < DIM; src_d++ ) {
            for ( int src_p = 0; src_p < 2; src_p++ ) {
               if (d_data[src_d][src_p].isInitialized()) {
                  dst_array.copy( d_data[src_d][src_p], box_list, src_offset );
               }
	    }
	 }
      }
   }
}


/*
*************************************************************************
*                                                                       *
* Print routines for outernode centered arrays.                         *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    std::ostream& os, 
                                    int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    int depth, 
                                    std::ostream& os, 
                                    int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array axis = " << axis << std::endl;
      for (int side = 0; side < 2; side++) {
         os << "Side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisSide(axis, side, box, depth, os, prec);
      }
   }
}

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::printAxisSide(int face_normal, 
                                            int side, 
                                            const hier::Box<DIM>& box, 
                                            std::ostream& os, 
                                            int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(face_normal, side, box, d, os, prec);
   }
}


template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::printAxisSide(int face_normal, 
                                            int side, 
                                            const hier::Box<DIM>& box, 
                                            int depth, 
                                            std::ostream& os, 
                                            int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   const hier::Box<DIM> nodebox = NodeGeometry<DIM>::toNodeBox(box);
   const hier::Box<DIM> region = nodebox * d_data[face_normal][side].getBox();
   os.precision(prec);
   for (typename hier::Box<DIM>::Iterator i(region); i; i++) {
      os << "array" << i() << " = " 
         << d_data[face_normal][side](i(),depth) << std::endl << std::flush;
   }
}

/*
*************************************************************************
*                                                                       *
* Checks that class version and restart file version are equal.         *
* If so, reads in d_depth from the database.                            *
* Then has each item in d_data read in its data from the database.      *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_OUTERNODEDATA_VERSION");
   if (ver != PDAT_OUTERNODEDATA_VERSION) {
      TBOX_ERROR("OuternodeData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
     sprintf(array_name, "d_data%d_1", i);
     if (database->keyExists(array_name)) {
        array_database = database->getDatabase(array_name);
        (d_data[i][0]).getFromDatabase(array_database);
     }

     sprintf(array_name, "d_data%d_2", i);
     if (database->keyExists(array_name)) {
        array_database = database->getDatabase(array_name);
        (d_data[i][1]).getFromDatabase(array_database);
     }
   }
}

/*
*************************************************************************
*                                                                       *
* Writes out class version number, d_depth to the database.             *
* Then has each item in d_data write out its data to the database.      *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuternodeData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_OUTERNODEDATA_VERSION",
                         PDAT_OUTERNODEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      if (d_data[i][0].isInitialized()) {
         sprintf(array_name, "d_data%d_1", i);
         array_database = database->putDatabase(array_name);
         (d_data[i][0]).putToDatabase(array_database);
      }
      if (d_data[i][1].isInitialized()) {
         sprintf(array_name, "d_data%d_2", i);
         array_database = database->putDatabase(array_name);
         (d_data[i][1]).putToDatabase(array_database);
      }
   }
}

}
}

#endif

