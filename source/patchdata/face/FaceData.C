//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/face/FaceData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated face centered patch data type
//

#ifndef included_pdat_FaceData_C
#define included_pdat_FaceData_C

#include "FaceData.h"

#include "Box.h"
#include "BoxList.h"
#include "FaceGeometry.h"
#include "FaceOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include <stdio.h>

#define PDAT_FACEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "FaceData.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for face data objects.  The constructor	*
* simply initializes data variables and sets up the array data.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
FaceData<DIM,TYPE>::FaceData(const hier::Box<DIM>& box,
                             int depth,
                             const hier::IntVector<DIM>& ghosts,
                             tbox::Pointer<tbox::Arena> pool)
:  hier::PatchData<DIM>(box, ghosts),
   d_depth(depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face = FaceGeometry<DIM>::toFaceBox(this -> getGhostBox(), d);
      d_data[d].initializeArray(face, depth, pool);
   }
}

template<int DIM, class TYPE>
FaceData<DIM,TYPE>::~FaceData()
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
FaceData<DIM,TYPE>::FaceData(const FaceData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::operator=(const FaceData<DIM,TYPE>& foo)
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const FaceData<DIM,TYPE> *t_src = 
      dynamic_cast<const FaceData<DIM,TYPE> *>(&src);

   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> box = d_data[d].getBox() * t_src->d_data[d].getBox();
         if (!box.empty()) {
            d_data[d].copy(t_src->d_data[d], box);
         }
      }
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   FaceData<DIM,TYPE> *t_dst =
      dynamic_cast<FaceData<DIM,TYPE> *>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> box = d_data[d].getBox() * t_dst->d_data[d].getBox();
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                              const hier::BoxOverlap<DIM>& overlap)
{
   const FaceData<DIM,TYPE> *t_src =
      dynamic_cast<const FaceData<DIM,TYPE> *>(&src);
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < DIM; d++) {
         hier::IntVector<DIM> face_offset(src_offset);
         if (d > 0) {
            for (int i = 0; i < DIM; i++) {
               face_offset(i) = src_offset((d+i)%DIM);
            }
         }
         const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList(d);
         d_data[d].copy(t_src->d_data[d], box_list, face_offset);
      }
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                               const hier::BoxOverlap<DIM>& overlap) const
{
   FaceData<DIM,TYPE> *t_dst = 
      dynamic_cast<FaceData<DIM,TYPE> *>(&dst);

   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::copyDepth(int dst_depth,
			           const FaceData<DIM,TYPE>& src,
				   int src_depth)
{
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> box = d_data[d].getBox() * src.d_data[d].getBox();
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

template<int DIM, class TYPE>
bool FaceData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int FaceData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < DIM; d++) {
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }
      size += d_data[d].getDataStreamSize(t_overlap->getDestinationBoxList(d),
                                          face_offset);
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
void FaceData<DIM,TYPE>::packStream(tbox::AbstractStream& stream,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }
      const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberItems() > 0) {
         d_data[d].packStream(stream, boxes, face_offset);
      }
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::unpackStream(tbox::AbstractStream& stream,
                                      const hier::BoxOverlap<DIM>& overlap)
{
   const FaceOverlap<DIM> *t_overlap =
      dynamic_cast<const FaceOverlap<DIM> *>(&overlap);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      hier::IntVector<DIM> face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < DIM; i++) {
            face_offset(i) = offset((d+i)%DIM);
         }
      }
      const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList(d);
      if (boxes.getNumberItems() > 0) {
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

template<int DIM, class TYPE>
size_t FaceData<DIM,TYPE>::getSizeOfData(const hier::Box<DIM>& box,
                                         int depth,
                                         const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   size_t size = 0;
   const hier::Box<DIM> ghost_box = hier::Box<DIM>::grow(box, ghosts);
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = FaceGeometry<DIM>::toFaceBox(ghost_box, d);
      size += ArrayData<DIM,TYPE>::getSizeOfData(face_box, depth);
   }
   return(size);
}

/*
*************************************************************************
*									*
* Fill the face centered box with the given value.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i].fill(t, d);
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::fill(const TYPE& t,
                              const hier::Box<DIM>& box,
                              int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i].fill(t, FaceGeometry<DIM>::toFaceBox(box, i), d);
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::fillAll(const TYPE& t)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i].fillAll(t);
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::fillAll(const TYPE& t, 
                                 const hier::Box<DIM>& box)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i].fillAll(t, FaceGeometry<DIM>::toFaceBox(box, i));
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                               std::ostream& os, 
                               int prec) const
{
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::print(const hier::Box<DIM>& box,
                               int d,
                               std::ostream& os,
                               int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int axis = 0; axis < DIM; axis++) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::printAxis(int face_normal,
                                   const hier::Box<DIM>& box,
                                   std::ostream& os,
                                   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxis(face_normal, box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void FaceData<DIM,TYPE>::printAxis(int face_normal, 
                                   const hier::Box<DIM>& box,
                                   int depth, 
                                   std::ostream& os,
                                   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   os.precision(prec);
   for (FaceIterator<DIM> i(box, face_normal); i; i++) {
      os << "array" << i() << " = " 
         << d_data[face_normal](i(),depth) << std::endl << std::flush;
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif
  
   int ver = database->getInteger("PDAT_FACEDATA_VERSION");
   if (ver != PDAT_FACEDATA_VERSION) {
      TBOX_ERROR("FaceData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   char array_name[8];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      sprintf(array_name, "d_data%d", i);
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

template<int DIM, class TYPE>
void FaceData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_FACEDATA_VERSION", PDAT_FACEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   char array_name[8];
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      sprintf(array_name, "d_data%d", i);
      array_database = database->putDatabase(array_name);
      (d_data[i]).putToDatabase(array_database);
   }
}



}
}

#endif
