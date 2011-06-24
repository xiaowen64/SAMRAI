/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Templated face centered patch data type 
 *
 ************************************************************************/

#ifndef included_pdat_FaceData_C
#define included_pdat_FaceData_C

#include "SAMRAI/pdat/FaceData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/FaceData.I"
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE> const int FaceData<TYPE>::PDAT_FACEDATA_VERSION = 1;

/*
 *************************************************************************
 *									*
 * Constructor and destructor for face data objects.  The constructor	*
 * simply initializes data variables and sets up the array data.		*
 *									*
 *************************************************************************
 */

template<class TYPE>
FaceData<TYPE>::FaceData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::Box face = FaceGeometry::toFaceBox(this->getGhostBox(), d);
      d_data[d].initializeArray(face, depth);
   }
}

template<class TYPE>
FaceData<TYPE>::~FaceData()
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

template<class TYPE>
FaceData<TYPE>::FaceData(
   const FaceData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<class TYPE>
void FaceData<TYPE>::operator = (
   const FaceData<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *									*
 * Perform a fast copy between two face centered arrays where their	*
 * index spaces overlap.							*
 *									*
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const FaceData<TYPE>* t_src =
      dynamic_cast<const FaceData<TYPE> *>(&src);

   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      for (int d = 0; d < getDim().getValue(); d++) {
         const hier::Box box = d_data[d].getBox() * t_src->d_data[d].getBox();
         if (!box.empty()) {
            d_data[d].copy(t_src->d_data[d], box);
         }
      }
   }
}

template<class TYPE>
void FaceData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   FaceData<TYPE>* t_dst =
      dynamic_cast<FaceData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != NULL);

   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::Box box = d_data[d].getBox() * t_dst->d_data[d].getBox();
      if (!box.empty()) {
         t_dst->d_data[d].copy(d_data[d], box);
      }
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

template<class TYPE>
void FaceData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const FaceData<TYPE>* t_src =
      dynamic_cast<const FaceData<TYPE> *>(&src);
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      const hier::IntVector& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); d++) {
         hier::IntVector face_offset(src_offset);
         if (d > 0) {
            for (int i = 0; i < getDim().getValue(); i++) {
               face_offset(i) = src_offset((d + i) % getDim().getValue());
            }
         }
         const hier::BoxList& box_list = t_overlap->getDestinationBoxList(d);
         d_data[d].copy(t_src->d_data[d], box_list, face_offset);
      }
   }
}

template<class TYPE>
void FaceData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   FaceData<TYPE>* t_dst =
      dynamic_cast<FaceData<TYPE> *>(&dst);

   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      hier::IntVector face_offset(src_offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = src_offset((d + i) % getDim().getValue());
         }
      }
      const hier::BoxList& box_list = t_overlap->getDestinationBoxList(d);
      t_dst->d_data[d].copy(d_data[d], box_list, face_offset);
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

template<class TYPE>
void FaceData<TYPE>::copyDepth(
   int dst_depth,
   const FaceData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::Box box = d_data[d].getBox() * src.d_data[d].getBox();
      if (!box.empty()) {
         d_data[d].copyDepth(dst_depth, src.d_data[d], src_depth, box);
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

template<class TYPE>
bool FaceData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int FaceData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < getDim().getValue(); d++) {
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      size += d_data[d].getDataStreamSize(t_overlap->getDestinationBoxList(d),
            face_offset);
   }
   return size;
}

/*
 *************************************************************************
 *									*
 * Pack/unpack data into/out of the message streams using the index	*
 * space in the overlap descriptor.					*
 *									*
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      const hier::BoxList& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberOfItems() > 0) {
         d_data[d].packStream(stream, boxes, face_offset);
      }
   }
}

template<class TYPE>
void FaceData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      const hier::BoxList& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberOfItems() > 0) {
         d_data[d].unpackStream(stream, boxes, face_offset);
      }
   }
}

/*
 *************************************************************************
 *									*
 * Calculate the amount of memory space needed to represent the data	*
 * for a  face centered grid.						*
 *									*
 *************************************************************************
 */

template<class TYPE>
size_t FaceData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   for (int d = 0; d < box.getDim().getValue(); d++) {
      const hier::Box face_box = FaceGeometry::toFaceBox(ghost_box, d);
      size += ArrayData<TYPE>::getSizeOfData(face_box, depth);
   }
   return size;
}

/*
 *************************************************************************
 *									*
 * Fill the face centered box with the given value.			*
 *									*
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i].fill(t, d);
   }
}

template<class TYPE>
void FaceData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i].fill(t, FaceGeometry::toFaceBox(box, i), d);
   }
}

template<class TYPE>
void FaceData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i].fillAll(t);
   }
}

template<class TYPE>
void FaceData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i].fillAll(t, FaceGeometry::toFaceBox(box, i));
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Print face centered data.  Note:  makes call to specialized printAxis *
 * routine in FaceDataSpecialized.C                                      *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<class TYPE>
void FaceData<TYPE>::print(
   const hier::Box& box,
   int d,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<class TYPE>
void FaceData<TYPE>::printAxis(
   int face_normal,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxis(face_normal, box, d, os, prec);
   }
}

template<class TYPE>
void FaceData<TYPE>::printAxis(
   int face_normal,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   os.precision(prec);
   for (FaceIterator i(box, face_normal); i; i++) {
      os << "array" << i() << " = "
      << d_data[face_normal](i(), depth) << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Checks that class version and restart file version are equal.  If so, *
 * reads in the d_depth data member to the database.  Then tells		*
 * d_data to read itself in from the database.				*
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
   TBOX_ASSERT(!database.isNull());

   int ver = database->getInteger("PDAT_FACEDATA_VERSION");
   if (ver != PDAT_FACEDATA_VERSION) {
      TBOX_ERROR("FaceData<getDim()>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = database->getDatabase(array_name);
      (d_data[i]).getFromDatabase(array_database);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Write out the class version number, d_depth data member to the	*
 * database.  Then tells d_data to write itself to the database.		*
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void FaceData<TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
   TBOX_ASSERT(!database.isNull());

   database->putInteger("PDAT_FACEDATA_VERSION", PDAT_FACEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = database->putDatabase(array_name);
      (d_data[i]).putToDatabase(array_database);
   }
}

}
}

#endif
