//
// File:	OuteredgeData.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 694 $
// Modified:	$Date: 2005-10-31 10:07:59 -0800 (Mon, 31 Oct 2005) $
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
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

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
   const int depth,
   tbox::Pointer<tbox::Arena> pool)
:  hier::PatchData<DIM>(box, hier::IntVector<DIM>(0)),
   d_depth(depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   // a = axis
   for (int a = 0; a < DIM; a++) {

      const hier::Box<DIM>& ghosts = this -> getGhostBox();
      hier::Box<DIM> edgebox = EdgeGeometry<DIM>::toEdgeBox(ghosts, a);

      // f = face normal
      for (int f = 0; f < DIM; f++) {

         // Data NULL when f == a
         if (f != a) {
            hier::Box<DIM> boxlo = edgebox;
            boxlo.upper(f) = edgebox.lower(f);
            hier::Box<DIM> boxup = edgebox;
            boxup.lower(f) = edgebox.upper(f);

            OuteredgeGeometry<DIM>::trimBoxes(boxlo, boxup, a, f);

            if (!boxlo.empty()) {
               d_data[a][f][0].initializeArray(boxlo, depth, pool);
            }
            if (!boxup.empty()) {
               d_data[a][f][1].initializeArray(boxup, depth, pool);
            }
            
         }
      }

   }
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
   }
   else if ( t_oedge_src != NULL ) {
      copyFromOuteredge( *t_oedge_src );
   }
   else {
      TBOX_ERROR("OuteredgeData can copy only with EdgeData<DIM> of the same"
		 <<"primitive data type.");
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
   }
   else if ( t_oedge_dst != NULL ) {
      copyToOuteredge( *t_oedge_dst );
   }
   else {
      TBOX_ERROR("OuteredgeData can copy only with EdgeData<DIM> of the same"
		 <<"primitive data type.");
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
   assert(t_overlap != NULL);
#endif

   const EdgeData<DIM,TYPE> *t_edge_src = 
      dynamic_cast<const EdgeData<DIM,TYPE> *>(&src);
   const OuteredgeData<DIM,TYPE> *t_oedge_src = 
      dynamic_cast<const OuteredgeData<DIM,TYPE> *>(&src);

   if ( t_edge_src != NULL ) {
      copyFromEdge( *t_edge_src, *t_overlap );
   }
   else if ( t_oedge_src != NULL ) {
      copyFromOuteredge( *t_oedge_src, *t_overlap );
   }
   else {
      TBOX_ERROR("OuteredgeData can copy only with EdgeData<DIM> of the same"
		 <<"primitive data type.");
   }

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif

   EdgeData<DIM,TYPE> *t_edge_dst = 
      dynamic_cast<EdgeData<DIM,TYPE> *>(&dst);
   OuteredgeData<DIM,TYPE> *t_oedge_dst = 
      dynamic_cast<OuteredgeData<DIM,TYPE> *>(&dst);

   if ( t_edge_dst != NULL ) {
      copyToEdge( *t_edge_dst, *t_overlap );
   }
   else if ( t_oedge_dst != NULL ) {
      copyToOuteredge( *t_oedge_dst, *t_overlap );
   }
   else {
      TBOX_ERROR("OuteredgeData can copy only with EdgeData<DIM> of the same"
		 <<"primitive data type.");
   }

}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy between two arrays at the                         *
* specified depths, where their	index spaces overlap.                   *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyDepth(int dst_depth,
                                          const EdgeData<DIM,TYPE>& src,
                                          int src_depth)
{
   // a = axis
   for (int a = 0; a < DIM; a++ ) {

      const ArrayData<DIM,TYPE> &edge_array = src.getArrayData(a);         

      // f = face normal
      for (int f = 0; f < DIM; f++) {

         // s = lower/upper 
         for ( int s = 0; s < 2; s++ ) {
            ArrayData<DIM,TYPE> &oedge_array = d_data[a][f][s];
            oedge_array.copyDepth(dst_depth,
                                  edge_array,
                                  src_depth,
                                  oedge_array.getBox());
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
   assert(t_overlap != NULL);
#endif
   int size = 0;

   // a = axis
   for (int a = 0; a < DIM; a++) {
      const hier::BoxList<DIM>& boxlist = 
         t_overlap->getDestinationBoxList(a);
      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();

      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            size += d_data[a][f][0].getDataStreamSize(boxlist, src_offset);
         }
         if (d_data[a][f][1].isInitialized()) {
            size += d_data[a][f][1].getDataStreamSize(boxlist, src_offset);
         }
      }
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
void OuteredgeData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif

   // a = axis
   for (int a = 0; a < DIM; a++) {
      const hier::BoxList<DIM>& dst_boxes    = t_overlap->getDestinationBoxList(a);
      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
      hier::IntVector<DIM> edge_offset(src_offset);
      
      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes); 
           dst_box; dst_box++) {
         const hier::Box<DIM> src_box = hier::Box<DIM>::shift(dst_box(), 
                                                              -edge_offset);
         // f = face normal
         for (int f = 0; f < DIM; f++) {
            
            // s = lower/upper 
            for (int s = 0; s < 2; s++) {
               const hier::Box<DIM> intersect = src_box * d_data[a][f][s].getBox();
               if (!intersect.empty()) {
                  d_data[a][f][s].packStream(stream, 
                                             hier::Box<DIM>::shift(intersect, 
                                                                   edge_offset),
                                             edge_offset);
               }
            }

         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const EdgeOverlap<DIM> *t_overlap =
      dynamic_cast<const EdgeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif

   // a = axis
   for (int a = 0; a < DIM; a++) {
      const hier::BoxList<DIM>& dst_boxes = 
         t_overlap->getDestinationBoxList(a);
      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
      for (typename hier::BoxList<DIM>::Iterator dst_box(dst_boxes); 
           dst_box; dst_box++) {

         // f = face normal
         for (int f = 0; f < DIM; f++) {
            for (int s = 0; s < 2; s++) {
               const hier::Box<DIM> intersect = 
                  dst_box() * d_data[a][f][s].getBox();
               if (!intersect.empty()) {
                  d_data[a][f][s].unpackStream(stream,intersect,src_offset);
               }
            }

         }
      }
   }
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
   const hier::Box<DIM>& box, const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   size_t size = 0;

   // a = axis
   for (int a = 0; a < DIM; a++) {
      hier::Box<DIM> loc0 = EdgeGeometry<DIM>::toEdgeBox(box,a);
      hier::Box<DIM> loc1 = EdgeGeometry<DIM>::toEdgeBox(box,a);

      // f = face normal
      for (int f = 0; f < DIM; f++) {

         if (f != a) {
            loc0.upper(f) = box.lower(f);
            loc1.lower(f) = box.upper(f);

            OuteredgeGeometry<DIM>::trimBoxes(loc0, loc1, a, f);
            
            size += ArrayData<DIM,TYPE>::getSizeOfData(loc0, depth)
               + ArrayData<DIM,TYPE>::getSizeOfData(loc1, depth);
         }
      }
   }
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
hier::Box<DIM> OuteredgeData<DIM,TYPE>::getDataBox( int axis, 
                                                    int face_nrml, 
                                                    int s )
{
   if ( axis < 0 || axis >=DIM || 
        face_nrml < 0 || face_nrml >=DIM || 
        s < 0 || s > 1 ) {
      TBOX_ERROR("Bad values for axis and/or face_nrml and/or s in\n"
                 "OuteredgeData<DIM>::getDataBox().\n");
   }

   /*
    * We start with the full box and chop it down to the databox
    * corresponding to the given dimension and s.
    */
   hier::Box<DIM> databox = EdgeGeometry<DIM>::toEdgeBox(this->getBox(),axis);
   
   // Data NULL when axis == face_nrml
   if (axis == face_nrml) {

      databox.setEmpty();

   } else {

      hier::Box<DIM> boxlo = databox;
      boxlo.upper(face_nrml) = databox.lower(face_nrml);
      hier::Box<DIM> boxup = databox;
      boxup.lower(face_nrml) = databox.upper(face_nrml);

      OuteredgeGeometry<DIM>::trimBoxes(boxlo, boxup, axis, face_nrml);
      
      if ( s == 0 ) {
         databox = boxlo;
      }
      else { // s == 1
         databox = boxup;
      }
   }
   return databox;
}

/*
*************************************************************************
*                                                                       *
* Fill the outeredge centered box with the given value.                 *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fill(const TYPE& t, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif

   // a = axis
   for (int a = 0; a < DIM; a++) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][0].fill(t, d);
         }
         if (d_data[a][f][1].isInitialized()) {
            d_data[a][f][1].fill(t, d);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fill(const TYPE& t,
                                const hier::Box<DIM>& box,
                                const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   // a = axis
   for (int a = 0; a < DIM; a++) {
      hier::Box<DIM> databox = EdgeGeometry<DIM>::toEdgeBox(box, a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][0].fill(t, databox, d);
         }
         if (d_data[a][f][1].isInitialized()) {
            d_data[a][f][1].fill(t, databox, d);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fillAll(const TYPE& t)
{
   // a = axis
   for (int a = 0; a < DIM; a++) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][0].fillAll(t);
         }
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][1].fillAll(t);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   // a = axis
   for (int a = 0; a < DIM; a++) {
      hier::Box<DIM> databox = EdgeGeometry<DIM>::toEdgeBox(box, a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][0].fillAll(t, databox);
         }
         if (d_data[a][f][1].isInitialized()) {
            d_data[a][f][1].fillAll(t, databox);
         }
      }
   }
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
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      const ArrayData<DIM,TYPE> &edge_array = src.getArrayData(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {
            if (d_data[a][f][s].isInitialized()) {
               ArrayData<DIM,TYPE> &oedge_array = d_data[a][f][s];
               oedge_array.copy( edge_array, oedge_array.getBox() );
            }
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToEdge(EdgeData<DIM,TYPE>& dst) const
{
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      ArrayData<DIM,TYPE> &edge_array = dst.getArrayData(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {
            if (d_data[a][f][s].isInitialized()) {
               edge_array.copy(d_data[a][f][s], 
                               d_data[a][f][s].getBox());
            }
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
void OuteredgeData<DIM,TYPE>::copyFromEdge(const EdgeData<DIM,TYPE>& src,
					     const EdgeOverlap<DIM>& overlap)
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            d_data[a][f][0].copy( src.getArrayData(a), 
                                  box_list, src_offset);
         }
         if (d_data[a][f][1].isInitialized()) {
            d_data[a][f][1].copy( src.getArrayData(a), 
                                  box_list, src_offset);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToEdge(EdgeData<DIM,TYPE>& dst,
					   const EdgeOverlap<DIM>& overlap) const
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (d_data[a][f][0].isInitialized()) {
            dst.getArrayData(a).copy(d_data[a][f][0], box_list, src_offset);
         }
         if (d_data[a][f][1].isInitialized()) {
            dst.getArrayData(a).copy(d_data[a][f][1], box_list, src_offset);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyFromOuteredge( 
   const OuteredgeData<DIM,TYPE> &src )
{
   // a = axis
   for (int a = 0; a < DIM; a++) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {
            if (src.d_data[a][f][s].isInitialized() && 
               d_data[a][f][s].isInitialized()) {
               const ArrayData<DIM,TYPE> &src_array = 
                  src.d_data[a][f][s];
               ArrayData<DIM,TYPE> &dst_array = 
                  d_data[a][f][s];
               dst_array.copy(src_array, 
                              dst_array.getBox() );
            }
         }
      }

   }

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyFromOuteredge( 
   const OuteredgeData<DIM,TYPE> &src,
   const EdgeOverlap<DIM> &overlap )
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      const hier::BoxList<DIM>& box_list =
         overlap.getDestinationBoxList(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {
            if (src.d_data[a][f][s].isInitialized() && 
               d_data[a][f][s].isInitialized()) {
               const ArrayData<DIM,TYPE> &src_array = 
                  src.d_data[a][f][s];
               d_data[a][f][s].copy(src_array, 
                                    box_list, src_offset);
            }
         }
      }

   } 

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToOuteredge( 
   OuteredgeData<DIM,TYPE> &dst ) const
{
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {
            if (dst.d_data[a][f][s].isInitialized() && 
                d_data[a][f][s].isInitialized()) {
               ArrayData<DIM,TYPE> &dst_array = dst.d_data[a][f][s];
               dst_array.copy( d_data[a][f][s],
                               d_data[a][f][s].getBox() );
            }
         }
      }   

   }

}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::copyToOuteredge( 
   OuteredgeData<DIM,TYPE> &dst,
   const EdgeOverlap<DIM> &overlap ) const
{
   const hier::IntVector<DIM>& src_offset = overlap.getSourceOffset();

   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      const hier::BoxList<DIM>& box_list = overlap.getDestinationBoxList(a);
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         for ( int s = 0; s < 2; s++ ) {

            if (dst.d_data[a][f][s].isInitialized() &&
                d_data[a][f][s].isInitialized()) {
               ArrayData<DIM,TYPE> &dst_array = dst.d_data[a][f][s];
               dst_array.copy( d_data[a][f][s], 
                               box_list, src_offset );
            }
         }
      }
   }

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
                                    ostream& os, 
                                    int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::print(
   const hier::Box<DIM>& box, const int d, ostream& os, int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   // a = axis
   for ( int a = 0; a < DIM; a++ ) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         os << "Array Axis Index<DIM> = " << a << "," << f << endl;
         for (int side = 0; side < 2; side++) {
            os << "Side Index<DIM> = " << ((side == 0) ? "lower" : "upper") 
               << endl;
            printAxisSide(a, f, side, box, d, os, prec);
         }
      }
   }
}

template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::printAxisSide(
   int axis, int face_nrml, int s, 
   const hier::Box<DIM>& box, ostream& os, int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM) && 
          (face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
#endif
   if (d_depth > 1) {
      for (int d = 0; d < d_depth; d++) {
         os << "Array Component Index<DIM> = " << d << endl;
         printAxisSide(axis, face_nrml, s, box, d, os, prec);
      }
   } else {
       printAxisSide(axis, face_nrml, s, box, 0, os, prec);
   }
}


template <int DIM, class TYPE>
void OuteredgeData<DIM,TYPE>::printAxisSide(
   int axis, int face_nrml, int s, 
   const hier::Box<DIM>& box, int d, ostream& os, int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
   assert((axis >= 0) && (axis < DIM));
   assert((face_nrml >= 0) && (face_nrml < DIM));
   assert((s == 0) || (s == 1));
#endif
   os.precision( ((prec < 0) ? 12 : prec) );
   if (axis != face_nrml) {
      const hier::Box<DIM> edgebox = EdgeGeometry<DIM>::toEdgeBox(box, axis);
      const hier::Box<DIM> region = edgebox * d_data[axis][face_nrml][s].getBox();
      os.precision( ((prec < 0) ? 12 : prec) );
      for (typename hier::Box<DIM>::Iterator ii(region); ii; ii++) {
         os << "array" << ii() << " = " << d_data[axis][face_nrml][s](ii(),d) << endl;
         os << flush;
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
   assert(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_OUTEREDGEDATA_VERSION");
   if (ver != PDAT_OUTEREDGEDATA_VERSION) {
      TBOX_ERROR("OuteredgeData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << endl);
   }

   d_depth = database->getInteger("d_depth");

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   // a = axis
   for (int a = 0; a < DIM; a++) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (f != a) {
            sprintf(array_name, "d_data%d_%d_1", a, f);
            array_database = database->getDatabase(array_name);
            (d_data[a][f][0]).getFromDatabase(array_database);
            
            sprintf(array_name, "d_data%d_%d_2", a, f);
            array_database = database->getDatabase(array_name);
            (d_data[a][f][1]).getFromDatabase(array_database);
         }
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
void OuteredgeData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif

   database->putInteger("PDAT_OUTEREDGEDATA_VERSION",
                         PDAT_OUTEREDGEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   // a = axis
   for (int a = 0; a < DIM; a++) {
      // f = face normal
      for (int f = 0; f < DIM; f++) {
         if (f != a) {
            sprintf(array_name, "d_data%d_%d_1", a, f);
            array_database = database->putDatabase(array_name);
            (d_data[a][f][0]).putToDatabase(array_database);
            
            sprintf(array_name, "d_data%d_%d_2", a, f);
            array_database = database->putDatabase(array_name);
            (d_data[a][f][1]).putToDatabase(array_database);
         }
      }
   }
}

}
}

#endif

