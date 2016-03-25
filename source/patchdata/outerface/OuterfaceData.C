//
// File:	OuterfaceData.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 179 $
// Modified:	$Date: 2005-01-20 14:50:51 -0800 (Thu, 20 Jan 2005) $
// Description:	Templated outerface centered patch data type
//

#ifndef included_pdat_OuterfaceData_C
#define included_pdat_OuterfaceData_C

#include "OuterfaceData.h"

#include "Box.h"
#include "BoxList.h"
#include "FaceData.h"
#include "FaceGeometry.h"
#include "FaceOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include <stdio.h>
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#define PDAT_OUTERFACEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "OuterfaceData.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for outerface data objects.  The		*
* constructor simply initializes data variables and sets up the		*
* array data.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OuterfaceData<DIM,TYPE>::OuterfaceData(
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

   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM>& ghosts = this -> getGhostBox();
      const hier::Box<DIM> facebox = FaceGeometry<DIM>::toFaceBox(ghosts, d);
      hier::Box<DIM> outerfacebox = facebox;
      outerfacebox.upper(0) = facebox.lower(0);
      d_data[d][0].initializeArray(outerfacebox, depth, pool);
      outerfacebox.lower(0) = facebox.upper(0);
      outerfacebox.upper(0) = facebox.upper(0);
      d_data[d][1].initializeArray(outerfacebox, depth, pool);
   }
}

template<int DIM, class TYPE>
OuterfaceData<DIM,TYPE>::~OuterfaceData()
{
}

/*
*************************************************************************
*									*
* The following are private and cannot be used, but they are defined	*
* here for compilers that require that every template declaration have	*
* a definition (a stupid requirement, if you ask me).			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OuterfaceData<DIM,TYPE>::OuterfaceData(
   const OuterfaceData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::operator=(const OuterfaceData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Perform a fast copy between an outerface patch data type (source) and	*
* a face patch data type (destination) where the index spaces overlap.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
    const FaceData<DIM,TYPE>* const t_src =
      dynamic_cast<const FaceData<DIM,TYPE> *>(&src);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_src != NULL);
#endif

   for (int axis = 0; axis < DIM; axis++ ) {
      const ArrayData<DIM,TYPE> &face_array = t_src->getArrayData(axis);
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE> &oface_array = d_data[axis][loc];
         oface_array.copy( face_array, oface_array.getBox() );
      }
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   FaceData<DIM,TYPE> *t_dst =
      dynamic_cast<FaceData<DIM,TYPE> *>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_dst != NULL);
#endif
   for (int d = 0; d < DIM; d++) {
      t_dst->getArrayData(d).copy(d_data[d][0], d_data[d][0].getBox());
      t_dst->getArrayData(d).copy(d_data[d][1], d_data[d][1].getBox());
   }
}

/*
*************************************************************************
*									*
* Copy data from the source into the destination according to the	*
* overlap descriptor.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                     const hier::BoxOverlap<DIM>& overlap)
{

   NULL_USE(src);
   NULL_USE(overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!"Copy with outerface as destination is not defined yet...");
#endif
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                      const hier::BoxOverlap<DIM>& overlap) const
{
   FaceData<DIM,TYPE> *t_dst =
      dynamic_cast<FaceData<DIM,TYPE> *>(&dst);
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_dst != NULL);
   assert(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      hier::IntVector<DIM> face_offset(src_offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = src_offset((d+i)%DIM);
         }
      }
      const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList(d);
      t_dst->getArrayData(d).copy(d_data[d][0], box_list, face_offset);
      t_dst->getArrayData(d).copy(d_data[d][1], box_list, face_offset);
   }
}

/*
*************************************************************************
*									*
* Perform a fast copy between two arrays at the                         *
* specified depths, where their	index spaces overlap.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::copyDepth(int dst_depth,
					  const FaceData<DIM,TYPE>& src,
					  int src_depth)
{
   for (int axis = 0; axis < DIM; axis++ ) {
      const ArrayData<DIM,TYPE> &face_array = src.getArrayData(axis);
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE> &oface_array = d_data[axis][loc];
	 oface_array.copyDepth(dst_depth,
			       face_array,
			       src_depth,
			       oface_array.getBox());
      }
   }
}

/*
*************************************************************************
*									*
* Calculate the buffer space needed to pack/unpack messages on the box	*
* region using the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool OuterfaceData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int OuterfaceData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxlist = t_overlap->getDestinationBoxList(d);
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }
      size += d_data[d][0].getDataStreamSize(boxlist, face_offset);
      size += d_data[d][1].getDataStreamSize(boxlist, face_offset);
   }
   return(size);
}

/*
*************************************************************************
*									*
* Pack/unpack data into/out of the message streams using the index	*
* space in the overlap descriptor.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const FaceOverlap<DIM> *t_overlap = 
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes    = t_overlap->getDestinationBoxList(d);
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }

      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         const hier::Box<DIM> src_box = hier::Box<DIM>::shift(b(), -face_offset);
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect = src_box * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].packStream(stream, 
                                       hier::Box<DIM>::shift(intersect, face_offset), 
                                       face_offset);
            }
         }
      }

   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes    = t_overlap->getDestinationBoxList(d);
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }

      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect = b() * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStream(stream, intersect, face_offset);
            }
         }
      }
   }
}

/*
*************************************************************************
*									*
* Calculate the amount of memory space needed to represent the data	*
* for a  outerface centered grid.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t OuterfaceData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box, const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   size_t size = 0;
   for (int d = 0; d < DIM; d++) {
      hier::Box<DIM> lower = FaceGeometry<DIM>::toFaceBox(box, d);
      hier::Box<DIM> upper = FaceGeometry<DIM>::toFaceBox(box, d);
      lower.upper(d) = box.lower(d);
      upper.lower(d) = box.upper(d);
      size += ArrayData<DIM,TYPE>::getSizeOfData(lower, depth);
      size += ArrayData<DIM,TYPE>::getSizeOfData(upper, depth);
   }
   return(size);
}

/*
*************************************************************************
*									*
* Fill the outerface centered box with the given value.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::fill(const TYPE& t, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fill(t, d);
      d_data[i][1].fill(t, d);
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::fill(const TYPE& t,
                                const hier::Box<DIM>& box,
                                const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fill(t, FaceGeometry<DIM>::toFaceBox(box, i), d);
      d_data[i][1].fill(t, FaceGeometry<DIM>::toFaceBox(box, i), d);
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::fillAll(const TYPE& t)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fillAll(t);
      d_data[i][1].fillAll(t);
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fillAll(t, FaceGeometry<DIM>::toFaceBox(box, i));
      d_data[i][1].fillAll(t, FaceGeometry<DIM>::toFaceBox(box, i));
   }
}

/*
*************************************************************************
*                                                                       *
* Print routines for outerface centered arrays.				*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::print(const hier::Box<DIM>& box, ostream& os, int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::print(
   const hier::Box<DIM>& box, const int d, ostream& os, int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array Axis hier::Index = " << axis << endl;
      for (int face = 0; face < 2; face++) {
         os << "Face hier::Index = " << ((face == 0) ? "lower" : "upper") << endl;
         printAxisFace(axis, face, box, d, os, prec);
      }
   }
}

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::printAxisFace(
   const int axis, const int face, const hier::Box<DIM>& box, ostream& os, int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
#endif
   if (d_depth > 1) {
      for (int d = 0; d < d_depth; d++) {
         os << "Array Component hier::Index = " << d << endl;
         printAxisFace(axis, face, box, d, os, prec);
      }
   } else {
       printAxisFace(axis, face, box, 0, os, prec);
   }
}

/*
*************************************************************************
*                                                                       *
* Checks that class version and restart file version are equal.  If so,	*
* reads in d_depth from the database.  Then has each item in d_data	*
* read in its data from the database.					*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_OUTERFACEDATA_VERSION");
   if (ver != PDAT_OUTERFACEDATA_VERSION) {
      TBOX_ERROR("OuterfaceData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << endl);
   }

   d_depth = database->getInteger("d_depth");

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
     sprintf(array_name, "d_data%d_1", i);
     array_database = database->getDatabase(array_name);
     (d_data[i][0]).getFromDatabase(array_database);

     sprintf(array_name, "d_data%d_2", i);
     array_database = database->getDatabase(array_name);
     (d_data[i][1]).getFromDatabase(array_database);
   }
}

/*
*************************************************************************
*                                                                       *
* Writes out class version number, d_depth to the database.		*
* Then has each item in d_data write out its data to the database.	*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif

   database->putInteger("PDAT_OUTERFACEDATA_VERSION",
                         PDAT_OUTERFACEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   char array_name[16];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
     sprintf(array_name, "d_data%d_1", i);
     array_database = database->putDatabase(array_name);
     (d_data[i][0]).putToDatabase(array_database);

     sprintf(array_name, "d_data%d_2", i);
     array_database = database->putDatabase(array_name);
     (d_data[i][1]).putToDatabase(array_database);
   }
}

/*
 * Print bool data
 */
template <int DIM, class TYPE>
void OuterfaceData<DIM,TYPE>::printAxisFace(
   const int axis, const int face,
   const hier::Box<DIM>& box, const int d, ostream& os,
   int prec) const
{
   NULL_USE(prec);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
#endif
   const hier::Box<DIM> facebox = FaceGeometry<DIM>::toFaceBox(box, axis);
   const hier::Box<DIM> region = facebox * d_data[axis][face].getBox();
   os.precision( ((prec < 0) ? 12 : prec) );
   for (typename hier::Box<DIM>::Iterator i(region); i; i++) {
      os << "array" << i() << " = " << d_data[axis][face](i(),d) << endl;
      os << flush;
   }
}


}
}

#endif
