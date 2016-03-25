//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerside/OutersideData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2085 $
// Modified:	$LastChangedDate: 2008-03-28 08:59:54 -0700 (Fri, 28 Mar 2008) $
// Description:	Templated outerside centered patch data type
//

#ifndef included_pdat_OutersideData_C
#define included_pdat_OutersideData_C

#include "OutersideData.h"

#include "Box.h"
#include "BoxList.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include <stdio.h>

#define PDAT_OUTERSIDEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "OutersideData.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for outerside data objects.  The		*
* constructor simply initializes data variables and sets up the		*
* array data.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OutersideData<DIM,TYPE>::OutersideData(
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
      const hier::Box<DIM>& ghosts = this -> getGhostBox();
      const hier::Box<DIM> sidebox = SideGeometry<DIM>::toSideBox(ghosts, d);
      hier::Box<DIM> outersidebox = sidebox;
      outersidebox.upper(d) = sidebox.lower(d);
      d_data[d][0].initializeArray(outersidebox, depth, pool);
      outersidebox.lower(d) = sidebox.upper(d);
      outersidebox.upper(d) = sidebox.upper(d);
      d_data[d][1].initializeArray(outersidebox, depth, pool);
   }
}

template<int DIM, class TYPE>
OutersideData<DIM,TYPE>::~OutersideData()
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
OutersideData<DIM,TYPE>::OutersideData(
   const OutersideData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::operator=(const OutersideData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Perform a fast copy between an outerside patch data type (source) and	*
* a side patch data type (destination) where the index spaces overlap.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
    const SideData<DIM,TYPE>* const t_src =
      dynamic_cast<const SideData<DIM,TYPE> *>(&src);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_src != NULL);
#endif

   for (int axis = 0; axis < DIM; axis++ ) {
      const ArrayData<DIM,TYPE> &side_array = t_src->getArrayData(axis);
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE> &oside_array = d_data[axis][loc];
         oside_array.copy( side_array, oside_array.getBox() );
      }
   }

}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   SideData<DIM,TYPE> *t_dst =
      dynamic_cast<SideData<DIM,TYPE> *>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
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
void OutersideData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                   const hier::BoxOverlap<DIM>& overlap)
{
   NULL_USE(src);
   NULL_USE(overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ERROR("Copy with outerside as destination is not defined yet...");
#endif
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   SideData<DIM,TYPE> *t_dst = 
      dynamic_cast<SideData<DIM,TYPE> *>(&dst);
   const SideOverlap<DIM> *t_overlap =
      dynamic_cast<const SideOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& box_list = t_overlap->getDestinationBoxList(d);
      t_dst->getArrayData(d).copy(d_data[d][0], box_list, src_offset);
      t_dst->getArrayData(d).copy(d_data[d][1], box_list, src_offset);
   }
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy from a side data object to this outerside data    *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::copyDepth(int dst_depth,
                                        const SideData<DIM,TYPE>& src,
                                        int src_depth)
{
   for (int axis = 0; axis < DIM; axis++ ) {
      const ArrayData<DIM,TYPE>& src_side_array = src.getArrayData(axis);
      for ( int loc = 0; loc < 2; loc++ ) {
         ArrayData<DIM,TYPE>& dst_oside_array = d_data[axis][loc];
	 dst_oside_array.copyDepth(dst_depth,
                                   src_side_array,
                                   src_depth,
                                   dst_oside_array.getBox());
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Perform a fast copy to a side data object from this outerside data    *
* object at the specified depths, where their index spaces overlap.     *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void OutersideData<DIM,TYPE>::copyDepth2(int dst_depth,
                                         SideData<DIM,TYPE>& dst,
                                         int src_depth) const
{
   for (int axis = 0; axis < DIM; axis++ ) {
      ArrayData<DIM,TYPE>& dst_side_array = dst.getArrayData(axis);
      for ( int loc = 0; loc < 2; loc++ ) {
         const ArrayData<DIM,TYPE>& src_oside_array = d_data[axis][loc];
         dst_side_array.copyDepth(dst_depth,
                                  src_oside_array,
                                  src_depth,
                                  src_oside_array.getBox());
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
bool OutersideData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int OutersideData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const SideOverlap<DIM> *t_overlap =
      dynamic_cast<const SideOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif

   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();  
 
   int size = 0;
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxlist = t_overlap->getDestinationBoxList(d);
      size += d_data[d][0].getDataStreamSize(boxlist, src_offset);
      size += d_data[d][1].getDataStreamSize(boxlist, src_offset);
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
void OutersideData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap) const
{
   const SideOverlap<DIM> *t_overlap =
      dynamic_cast<const SideOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes    = t_overlap->getDestinationBoxList(d);
      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         const hier::Box<DIM> src_box = hier::Box<DIM>::shift(b(), -src_offset);
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect = src_box * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].packStream(stream, 
                                       hier::Box<DIM>::shift(intersect, src_offset), 
                                       src_offset);
            }
         }
      }
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxOverlap<DIM>& overlap)
{
   const SideOverlap<DIM> *t_overlap =
      dynamic_cast<const SideOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   const hier::IntVector<DIM>& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < DIM; d++) {
      const hier::BoxList<DIM>& boxes = t_overlap->getDestinationBoxList(d);
      for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
         for (int f = 0; f < 2; f++) {
            const hier::Box<DIM> intersect = b() * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStream(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
*************************************************************************
*									*
* Calculate the amount of memory space needed to represent the data	*
* for a  outerside centered grid.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t OutersideData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   size_t size = 0;
   for (int d = 0; d < DIM; d++) {
      hier::Box<DIM> lower = SideGeometry<DIM>::toSideBox(box, d);
      hier::Box<DIM> upper = SideGeometry<DIM>::toSideBox(box, d);
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
* Fill the outerside centered box with the given value.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fill(t, d);
      d_data[i][1].fill(t, d);
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::fill(const TYPE& t,
                                   const hier::Box<DIM>& box,
                                   int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fill(t, SideGeometry<DIM>::toSideBox(box, i), d);
      d_data[i][1].fill(t, SideGeometry<DIM>::toSideBox(box, i), d);
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::fillAll(const TYPE& t)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fillAll(t);
      d_data[i][1].fillAll(t);
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::fillAll(const TYPE& t, 
                                      const hier::Box<DIM>& box)
{
   for (int i = 0; i < DIM; i++) {
      d_data[i][0].fillAll(t, SideGeometry<DIM>::toSideBox(box, i));
      d_data[i][1].fillAll(t, SideGeometry<DIM>::toSideBox(box, i));
   }
}

/*
*************************************************************************
*                                                                       *
* Print routines for outerside centered arrays.				*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    std::ostream& os, 
                                    int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                                    int depth, 
                                    std::ostream& os,
                                    int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   for (int side_normal = 0; side_normal < DIM; side_normal++) {
      os << "Array side normal  = " << side_normal << std::endl;
      for (int side = 0; side < 2; side++) {
         os << "side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisSide(side_normal, side, box, depth, os, prec);
      }
   }
}

template<int DIM, class TYPE>
void OutersideData<DIM,TYPE>::printAxisSide(int side_normal, 
                                            int side, 
                                            const hier::Box<DIM>& box, 
                                            std::ostream& os, 
                                            int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(side_normal, side, box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void OutersideData<DIM,TYPE>::printAxisSide(int side_normal, 
                                            int side,
                                            const hier::Box<DIM>& box,                                                               int depth,                     
                                            std::ostream& os,
                                            int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   const hier::Box<DIM> sidebox = 
      SideGeometry<DIM>::toSideBox(box, side_normal);
   const hier::Box<DIM> region = 
      sidebox * d_data[side_normal][side].getBox();
   os.precision(prec);
   for (typename hier::Box<DIM>::Iterator i(region); i; i++) {
      os << "array" << i() << " = " 
         << d_data[side_normal][side](i(),depth) << std::endl;
      os << std::flush;
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
void OutersideData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_OUTERSIDEDATA_VERSION");
   if (ver != PDAT_OUTERSIDEDATA_VERSION) {
      TBOX_ERROR("OutersideData<DIM>::getSpecializedFromDatabase error...\n"
		 << " : Restart file version different than class version" << std::endl);
   }
   
   d_depth = database->getInteger("d_depth");
   
   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i) + "_1";
      array_database = database->getDatabase(array_name);
      (d_data[i][0]).getFromDatabase(array_database);
      
      array_name = "d_data%d_" + tbox::Utilities::intToString(i) + "_2";
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
void OutersideData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_OUTERSIDEDATA_VERSION",
                         PDAT_OUTERSIDEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   for (int i = 0; i < DIM; i++) {
      std::string array_name = "d_data%d_" + tbox::Utilities::intToString(i) + "_1";
      array_database = database->putDatabase(array_name);
      (d_data[i][0]).putToDatabase(array_database);
      
      array_name = "d_data%d_" + tbox::Utilities::intToString(i) + "_2";
      array_database = database->putDatabase(array_name);
      (d_data[i][1]).putToDatabase(array_database);
   }
}


}
}

#endif
