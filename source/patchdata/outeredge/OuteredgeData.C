//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/patchdata/outeredge/OuteredgeData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Templated outeredge centered patch data type
//

#ifndef included_pdat_OuteredgeData_C
#define included_pdat_OuteredgeData_C


#include "OuteredgeData.h"

#include "Box.h"
#include "BoxList.h"
#include "EdgeData.h"
#include "EdgeGeometry.h"
#include "EdgeOverlap.h"
#include "OuteredgeGeometry.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include <stdio.h>

#define PDAT_OUTEREDGEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "OuteredgeData.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for outeredge data objects.  The           *
* constructor simply initializes data variables and sets up the         *
* array data.                                                           *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuteredgeData<DIM,TYPE>::OuteredgeData(
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

   for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                 hier::Box<DIM> oedge_data_box = 
                    OuteredgeGeometry<DIM>::toOuteredgeBox(this->getGhostBox(),
                                                           axis,
                                                           face_normal,
                                                           side);

                 if ( !oedge_data_box.empty() ) {
                    d_data[axis][face_normal][side].
                       initializeArray(oedge_data_box, depth, pool);
                 } 

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
OuteredgeData<DIM,TYPE>::~OuteredgeData()
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
OuteredgeData<DIM,TYPE>::OuteredgeData(
   const OuteredgeData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::operator=(const OuteredgeData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy between an outeredge patch data type (source) and *
* a edge patch data type (destination) where the index spaces overlap.  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const EdgeData<DIM,TYPE>* const t_edge_src =
      dynamic_cast<const EdgeData<DIM,TYPE> *>(&src);
   const OuteredgeData<DIM,TYPE>* const t_oedge_src =
      dynamic_cast<const OuteredgeData<DIM,TYPE> *>(&src);

   if ( t_edge_src != NULL ) {
      copyFromEdge( *t_edge_src );
   } else if ( t_oedge_src != NULL ) {
      copyFromOuteredge( *t_oedge_src );
   } else {
      TBOX_ERROR("OuteredgeData<DIM>::copy error!\n"
                 << "Can copy only from EdgeData<DIM,TYPE> or "
                 << "OuteredgeData<DIM,TYPE> of the same "
                 << "DIM and TYPE.");
   }

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   EdgeData<DIM,TYPE> *t_edge_dst =
      dynamic_cast<EdgeData<DIM,TYPE> *>(&dst);
   OuteredgeData<DIM,TYPE> *t_oedge_dst =
      dynamic_cast<OuteredgeData<DIM,TYPE> *>(&dst);

   if ( t_edge_dst != NULL ) {
      copyToEdge( *t_edge_dst );
   } else if ( t_oedge_dst != NULL ) {
      copyToOuteredge( *t_oedge_dst );
   } else {
      TBOX_ERROR("OuteredgeData<DIM>::copy2 error!\n"
                 << "Can copy only to EdgeData<DIM,TYPE> or "
                 << "OuteredgeData<DIM,TYPE> of the same "
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
void OuteredgeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                   const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const EdgeData<DIM,TYPE> *t_edge_src = 
      dynamic_cast<const EdgeData<DIM,TYPE> *>(&src);
   const OuteredgeData<DIM,TYPE> *t_oedge_src = 
      dynamic_cast<const OuteredgeData<DIM,TYPE> *>(&src);

   if ( t_edge_src != NULL ) {
      copyFromEdge( *t_edge_src, *t_overlap );
   } else if ( t_oedge_src != NULL ) {
      copyFromOuteredge( *t_oedge_src, *t_overlap );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy error!\n"
                 << "Can copy only from EdgeData<DIM,TYPE> or "
                 << "OuteredgeData<DIM,TYPE> of the same "
                 << "DIM and TYPE.");
   }

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   EdgeData<DIM,TYPE> *t_edge_dst = 
      dynamic_cast<EdgeData<DIM,TYPE> *>(&dst);
   OuteredgeData<DIM,TYPE> *t_oedge_dst = 
      dynamic_cast<OuteredgeData<DIM,TYPE> *>(&dst);

   if ( t_edge_dst != NULL ) {
      copyToEdge( *t_edge_dst, *t_overlap );
   } else if ( t_oedge_dst != NULL ) {
      copyToOuteredge( *t_oedge_dst, *t_overlap );
   } else {
      TBOX_ERROR("OuternodeData<DIM>::copy2 error!\n"
                 << "Can copy only to EdgeData<DIM,TYPE> or "
                 << "OuteredgeData<DIM,TYPE> of the same "
                 << "DIM and TYPE.");
   }

}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy from an edge data object to this outeredge data   *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyDepth(int dst_depth,
                                        const EdgeData<DIM,TYPE>& src,
                                        int src_depth)
{

  for (int axis = 0; axis < DIM; ++axis) {

     const ArrayData<DIM,TYPE>& src_edge_array = src.getArrayData(axis);

     for (int face_normal = 0; face_normal < DIM; ++face_normal) {

        if ( face_normal != axis ) {

           for (int side = 0; side < 2; ++side) {

              ArrayData<DIM,TYPE>& dst_oedge_array = 
                 d_data[axis][face_normal][side];

              dst_oedge_array.copyDepth(dst_depth,
                                        src_edge_array,
                                        src_depth,
                                        dst_oedge_array.getBox());

           }  // iterate over lower/upper sides

        }  // data is undefined when axis == face_normal

     }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy to an edge data object from this outeredge data   *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyDepth2(int dst_depth,
                                         EdgeData<DIM,TYPE>& dst,
                                         int src_depth) const
{

   for (int axis = 0; axis < DIM; ++axis) {

      ArrayData<DIM,TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               const ArrayData<DIM,TYPE>& src_oedge_array =
                  d_data[axis][face_normal][side];

               dst_edge_array.copyDepth(dst_depth,
                                        src_oedge_array,
                                        src_depth,
                                        src_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
*************************************************************************
*                                                                       *
* Add source data to the destination according to overlap descriptor.   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::sum(
   const hier::PatchData<DIM>& src,
   const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const OuteredgeData<DIM,TYPE> *t_oedge_src =
      dynamic_cast<const OuteredgeData<DIM,TYPE> *>(&src);

   // NOTE:  We assume this operation is only needed to
   //        copy and add data to another outeredge data
   //        object.  If we ever need to provide this for edge
   //        data or other flavors of the copy operation, we
   //        should refactor the routine similar to the way
   //        the regular copy operations are implemented.
   if ( t_oedge_src == NULL ) {
      TBOX_ERROR("OuteredgeData<DIM>::sum error!\n"
                 << "Can copy and add only from OuteredgeData<DIM,TYPE> "
                 << "of the same DIM and TYPE.");
   } else {

      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

      for (int axis = 0; axis < DIM; ++axis) {

         const hier::BoxList<DIM>& box_list =
            t_overlap->getDestinationBoxList(axis);

         for (int src_face_normal = 0; src_face_normal < DIM; ++src_face_normal) {

            if ( src_face_normal != axis ) {

               for (int src_side = 0; src_side < 2; ++src_side) {

                  if (t_oedge_src->d_data[axis][src_face_normal][src_side].
                         isInitialized() ) {

                     const ArrayData<DIM,TYPE>& src_array =
                        t_oedge_src->d_data[axis][src_face_normal][src_side];

                     for (int dst_face_normal = 0; 
                          dst_face_normal < DIM; ++dst_face_normal) {

                        if ( dst_face_normal != axis ) {

                           for (int dst_side = 0; dst_side < 2; ++dst_side) { 
                      
                              if (d_data[axis][dst_face_normal][dst_side].
                                     isInitialized() ) {

                                 d_data[axis][dst_face_normal][dst_side].
                                    sum(src_array,
                                        box_list, 
                                        src_offset);

                              }  // if dst data array is initialized

                           }  // iterate over dst lower/upper sides

                        }  // dst data is undefined when axis == face_normal

                     }  // iterate over dst face normal directions

                  }  // if src data array is initialized

               }  // iterate over src lower/upper sides

            }  // src data is undefined when axis == face_normal

         }  // iterate over src face normal directions

      }  // iterate over axis directions

   } // else t_oedge_src != NULL

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
bool OuteredgeData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template <int DIM, class TYPE>
int OuteredgeData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   int size = 0;

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& boxlist =
         t_overlap->getDestinationBoxList(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               if ( d_data[axis][face_normal][side].isInitialized() ) { 

                  size += d_data[axis][face_normal][side].
                                getDataStreamSize(boxlist, src_offset); 

                }  // if data arrays is initialized

             }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

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
void OuteredgeData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& dst_boxes =
         t_overlap->getDestinationBoxList(axis);

      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes);
           dst_box; dst_box++) {

         const hier::Box<DIM> src_box = hier::Box<DIM>::shift(dst_box(),
                                                              -src_offset);

         for (int face_normal = 0; face_normal < DIM; ++face_normal) {

            if ( face_normal != axis ) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box<DIM> intersection = 
                     src_box * d_data[axis][face_normal][side].getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side].
                        packStream(stream, 
                                   hier::Box<DIM>::shift(intersection, 
                                                         src_offset),
                                   src_offset);

                  } // if intersection non-empty
 
               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& dst_boxes =
         t_overlap->getDestinationBoxList(axis);

      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes);
           dst_box; dst_box++) {

         for (int face_normal = 0; face_normal < DIM; ++face_normal) {

            if ( face_normal != axis ) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box<DIM> intersection =
                     dst_box() * d_data[axis][face_normal][side].getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side].unpackStream(stream, 
                                                                  intersection, 
                                                                  src_offset);

                  } // if intersection non-empty

               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

/*
*************************************************************************
*                                                                       *
* Unpack data from the message stream and add to this outeredge data    *
* object using the index space in the overlap descriptor.               *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::unpackStreamAndSum(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& dst_boxes =
         t_overlap->getDestinationBoxList(axis);

      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes);
           dst_box; dst_box++) {

         for (int face_normal = 0; face_normal < DIM; ++face_normal) {

            if ( face_normal != axis ) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box<DIM> intersection =
                     dst_box() * d_data[axis][face_normal][side].getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side].
                        unpackStreamAndSum(stream, 
                                           intersection, 
                                           src_offset);

                  } // if intersection non-empty

               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

/*
*************************************************************************
*                                                                       *
* Calculate the amount of memory space needed to represent the data     *
* for a  outeredge centered grid.                                       *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
size_t OuteredgeData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box, 
   int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   size_t size = 0;

   for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                 hier::Box<DIM> oedge_data_box =
                    OuteredgeGeometry<DIM>::toOuteredgeBox(box,
                                                           axis,
                                                           face_normal,
                                                           side);

                 size += 
                    ArrayData<DIM,TYPE>::getSizeOfData(oedge_data_box, depth); 

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

   return(size);
}

/*
*************************************************************************
*                                                                       *
* Compute the box of valid edge indices given values of                 *
* dimension and side designating the set of data indices.               *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
hier::Box<DIM> 
OuteredgeData<DIM,TYPE>::getDataBox(int axis, 
                                    int face_normal, 
                                    int side )
{
   return(
      OuteredgeGeometry<DIM>::toOuteredgeBox(this->getGhostBox(),
                                             axis,
                                             face_normal,
                                             side)
   ); 
}

/*
*************************************************************************
*                                                                       *
* Fill the outeredge centered box with the given value.                 *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fill(const TYPE& t, 
                                   int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif

   for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                if ( d_data[axis][face_normal][side].isInitialized()) {
                   d_data[axis][face_normal][side].fill(t, d);
                }

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fill(const TYPE& t,
                                   const hier::Box<DIM>& box,
                                   int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif

   for (int axis = 0; axis < DIM; ++axis) {

      hier::Box<DIM> databox = EdgeGeometry<DIM>::toEdgeBox(box, axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                if ( d_data[axis][face_normal][side].isInitialized()) {
                   d_data[axis][face_normal][side].fill(t, databox, d);
                }

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fillAll(const TYPE& t)
{

   for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                if ( d_data[axis][face_normal][side].isInitialized()) {
                   d_data[axis][face_normal][side].fillAll(t);
                }

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fillAll(const TYPE& t, 
                                      const hier::Box<DIM>& box)
{

   for (int axis = 0; axis < DIM; ++axis) {

      hier::Box<DIM> databox = EdgeGeometry<DIM>::toEdgeBox(box, axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

          if ( face_normal != axis ) {

             for (int side = 0; side < 2; ++side) {

                if ( d_data[axis][face_normal][side].isInitialized()) {
                   d_data[axis][face_normal][side].fillAll(t, databox);
                }

             }  // iterate over lower/upper sides

          }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy between an outeredge patch data type (source) and *
* a edge patch data type (destination) where the index spaces overlap.  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyFromEdge(const EdgeData<DIM,TYPE>& src)
{

   for (int axis = 0; axis < DIM; ++axis) {

      const ArrayData<DIM,TYPE>& src_edge_array = src.getArrayData(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {
 
            for (int side = 0; side < 2; ++side) {

               ArrayData<DIM,TYPE>& dst_oedge_array =
                  d_data[axis][face_normal][side];

               dst_oedge_array.copy(src_edge_array,
                                    dst_oedge_array.getBox());

            }  // iterate over lower/upper sides
 
         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions
 
   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToEdge(EdgeData<DIM,TYPE>& dst) const
{

   for (int axis = 0; axis < DIM; ++axis) {

      ArrayData<DIM,TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {
 
            for (int side = 0; side < 2; ++side) {

               const ArrayData<DIM,TYPE>& src_oedge_array =
                  d_data[axis][face_normal][side];

               dst_edge_array.copy(src_oedge_array,
                                   src_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

    }  // iterate over axis directions

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
void OuteredgeData<DIM,TYPE>::copyFromEdge(const EdgeData<DIM,TYPE>& src,
                                           const EdgeOverlap<DIM>& overlap)
{

   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList(axis);
      const ArrayData<DIM,TYPE>& src_edge_array = src.getArrayData(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               ArrayData<DIM,TYPE>& dst_oedge_array =
                  d_data[axis][face_normal][side];
 
               dst_oedge_array.copy(src_edge_array,
                                    box_list,
                                    src_offset);

            }  // iterate over lower/upper sides
 
         }  // data is undefined when axis == face_normal
 
      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToEdge(
   EdgeData<DIM,TYPE>& dst,
   const EdgeOverlap<DIM>& overlap) const
{

   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList(axis);
      ArrayData<DIM,TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               const ArrayData<DIM,TYPE>& src_oedge_array =
                  d_data[axis][face_normal][side];

               dst_edge_array.copy(src_oedge_array,
                                   box_list,
                                   src_offset);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

    }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyFromOuteredge( 
   const OuteredgeData<DIM,TYPE>& src )
{

   for (int axis = 0; axis < DIM; ++axis) {

      for (int src_face_normal = 0; src_face_normal < DIM; ++src_face_normal) {

         if ( src_face_normal != axis ) {

            for (int src_side = 0; src_side < 2; ++src_side) {

               const ArrayData<DIM,TYPE>& src_oedge_array =
                  src.d_data[axis][src_face_normal][src_side];

               for (int dst_face_normal = 0; dst_face_normal < DIM; ++dst_face_normal) {

                  if ( dst_face_normal != axis ) {

                     for (int dst_side = 0; dst_side < 2; ++dst_side) {

                        ArrayData<DIM,TYPE>& dst_oedge_array =
                           d_data[axis][dst_face_normal][dst_side];

                        dst_oedge_array.copy(src_oedge_array,
                                             dst_oedge_array.getBox());

                     }  // iterate over dst lower/upper sides

                  }  // dst data is undefined when axis == face_normal
  
               }  // iterate over dst face normal directions

            }  // iterate over src lower/upper sides

         }  // src data is undefined when axis == face_normal
  
      }  // iterate over src face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyFromOuteredge( 
   const OuteredgeData<DIM,TYPE>& src,
   const EdgeOverlap<DIM>& overlap )
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < DIM; ++axis) {

      const hier::BoxList<DIM>& box_list =
         overlap.getDestinationBoxList(axis);

      for (int src_face_normal = 0; src_face_normal < DIM; ++src_face_normal) {

         if ( src_face_normal != axis ) {

            for (int src_side = 0; src_side < 2; ++src_side) {

               const ArrayData<DIM,TYPE>& src_oedge_array =
                  src.d_data[axis][src_face_normal][src_side];

               for (int dst_face_normal = 0; dst_face_normal < DIM; ++dst_face_normal) {

                  if ( dst_face_normal != axis ) {

                     for (int dst_side = 0; dst_side < 2; ++dst_side) {

                        ArrayData<DIM,TYPE>& dst_oedge_array =
                           d_data[axis][dst_face_normal][dst_side];

                        dst_oedge_array.copy(src_oedge_array,
                                             box_list,
                                             src_offset);

                     }  // iterate over dst lower/upper sides

                  }  // dst data is undefined when axis == face_normal
  
               }  // iterate over dst face normal directions

            }  // iterate over src lower/upper sides

         }  // src data is undefined when axis == face_normal
  
      }  // iterate over src face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToOuteredge( 
   OuteredgeData<DIM,TYPE>& dst ) const
{
   dst.copyFromOuteredge( *this );
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToOuteredge( 
   OuteredgeData<DIM,TYPE>& dst,
   const EdgeOverlap<DIM>& overlap ) const
{
   dst.copyFromOuteredge( *this, overlap );
}


/*
*************************************************************************
*                                                                       *
* Print routines for outeredge centered arrays.                         *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    std::ostream& os, 
                                    int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    int depth, 
                                    std::ostream& os, 
                                    int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif

  for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         os << "Array axis, face normal = " 
            << axis << "," << face_normal << std::endl;

         for (int side = 0; side < 2; ++side) {

            os << "side  = " 
               << ((side == 0) ? "lower" : "upper") << std::endl;

            printAxisSide(axis, face_normal, side, 
                          box, depth, os, prec);

         }  // iterate over lower/upper sides

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::printAxisSide(
   int axis, 
   int face_normal, 
   int side, 
   const hier::Box<DIM>& box, 
   std::ostream& os, 
   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(axis, face_normal, side, 
                    box, d, os, prec);
   }

}


template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::printAxisSide(
   int axis, 
   int face_normal, 
   int side, 
   const hier::Box<DIM>& box, 
   int depth, 
   std::ostream& os, 
   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif

   NULL_USE(prec);

   if (axis == face_normal) {

      os << "array data undefined" << std::endl;

   } else {

      const hier::Box<DIM> edgebox = EdgeGeometry<DIM>::toEdgeBox(box, axis);
      const hier::Box<DIM> region = 
         edgebox * d_data[axis][face_normal][side].getBox();
      for (typename hier::Box<DIM>::Iterator ii(region); ii; ii++) {
         os << "array" << ii() << " = " 
            << d_data[axis][face_normal][side](ii(),depth) << std::endl;
         os << std::flush;
      }

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
void OuteredgeData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_OUTEREDGEDATA_VERSION");
   if (ver != PDAT_OUTEREDGEDATA_VERSION) {
      TBOX_ERROR("OuteredgeData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;

   for (int axis = 0; axis < DIM; ++axis) {
 
      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               sprintf(array_name, "d_data%d_%d_%d", axis, face_normal, side);
               array_database = database->getDatabase(array_name);
               (d_data[axis][face_normal][side]).getFromDatabase(array_database);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

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
void OuteredgeData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_OUTEREDGEDATA_VERSION",
                         PDAT_OUTEREDGEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;

   for (int axis = 0; axis < DIM; ++axis) {

      for (int face_normal = 0; face_normal < DIM; ++face_normal) {

         if ( face_normal != axis ) {

            for (int side = 0; side < 2; ++side) {

               sprintf(array_name, "d_data%d_%d_%d", axis, face_normal, side);
               array_database = database->putDatabase(array_name);
               (d_data[axis][face_normal][side]).putToDatabase(array_database);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

}
}

#endif

