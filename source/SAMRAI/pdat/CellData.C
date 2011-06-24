/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Templated cell centered patch data type 
 *
 ************************************************************************/

#ifndef included_pdat_CellData_C
#define included_pdat_CellData_C

#include "SAMRAI/pdat/CellData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/CellData.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE> const int CellData<TYPE>::PDAT_CELLDATA_VERSION = 1;

/*
 *************************************************************************
 *									*
 * Calculate the amount of memory space needed to represent the data	*
 * for a cell centered grid.						*
 *									*
 *************************************************************************
 */

template<class TYPE>
size_t CellData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);

   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   return ArrayData<TYPE>::getSizeOfData(ghost_box, depth);
}

/*
 *************************************************************************
 *									*
 * Constructor and destructor for cell data objects.  The constructor	*
 * simply initializes data variables and sets up the array data.		*
 *									*
 *************************************************************************
 */

template<class TYPE>
CellData<TYPE>::CellData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth),
   d_data()
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   d_data.initializeArray(this->getGhostBox(), depth);
}

template<class TYPE>
CellData<TYPE>::~CellData()
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
CellData<TYPE>::CellData(
   const CellData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth()),
   d_data()
{
   NULL_USE(foo);
}

template<class TYPE>
void CellData<TYPE>::operator = (
   const CellData<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *									*
 * Perform a fast copy between two cell centered arrays where their	*
 * index spaces overlap.							*
 *									*
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_data, src);

   const CellData<TYPE>* t_src =
      dynamic_cast<const CellData<TYPE> *>(&src);
   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      const hier::Box box = d_data.getBox() * t_src->d_data.getBox();
      if (!box.empty()) {
         d_data.copy(t_src->d_data, box);
      }
   }
}

template<class TYPE>
void CellData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_data, dst);

   CellData<TYPE>* t_dst = dynamic_cast<CellData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != NULL);

   const hier::Box box = d_data.getBox() * t_dst->d_data.getBox();
   if (!box.empty()) {
      t_dst->d_data.copy(d_data, box);
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
void CellData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   const CellData<TYPE>* t_src =
      dynamic_cast<const CellData<TYPE> *>(&src);

   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      d_data.copy(t_src->d_data,
         t_overlap->getDestinationBoxList(),
         t_overlap->getSourceOffset());
   }
}

template<class TYPE>
void CellData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   CellData<TYPE>* t_dst = dynamic_cast<CellData<TYPE> *>(&dst);
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   t_dst->d_data.copy(d_data,
      t_overlap->getDestinationBoxList(),
      t_overlap->getSourceOffset());
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
void CellData<TYPE>::copyDepth(
   int dst_depth,
   const CellData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_data, src);

   const hier::Box box = d_data.getBox() * src.d_data.getBox();
   if (!box.empty()) {
      d_data.copyDepth(dst_depth, src.d_data, src_depth, box);
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
bool CellData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int CellData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   return d_data.getDataStreamSize(t_overlap->getDestinationBoxList(),
      t_overlap->getSourceOffset());
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
void CellData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   d_data.packStream(stream, t_overlap->getDestinationBoxList(),
      t_overlap->getSourceOffset());
}

template<class TYPE>
void CellData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   d_data.unpackStream(stream, t_overlap->getDestinationBoxList(),
      t_overlap->getSourceOffset());
}

/*
 *************************************************************************
 *									*
 * Print cell centered data.  Note:  makes call to specialized print     *
 * routine in CellDataSpecialized.C                                     *
 *									*
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      print(box, d, os, prec);
   }
}

template<class TYPE>
void CellData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   os.precision(prec);
   for (CellIterator i(box); i; i++) {
      os << "array" << i() << " = "
         << d_data(i(), depth) << std::endl << std::flush;
      os << std::flush;
   }
}

/*
 *************************************************************************
 *									*
 * Checks that class version and restart file version are equal.  If so,	*
 * reads in the d_depth data member to the database.  Then tells		*
 * d_data to read itself in from the database.				*
 *									*
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{

   TBOX_ASSERT(!database.isNull());

   int ver = database->getInteger("PDAT_CELLDATA_VERSION");
   if (ver != PDAT_CELLDATA_VERSION) {
      TBOX_ERROR("CellData<DIM>::getSpecializedFromDatabase error...\n"
         << "Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   tbox::Pointer<tbox::Database> array_database;
   array_database = database->getDatabase("d_data");
   d_data.getFromDatabase(array_database);
}

/*
 *************************************************************************
 *									*
 * Write out the class version number, d_depth data member to the	*
 * database.  Then tells d_data to write itself to the database.		*
 *									*
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
   TBOX_ASSERT(!database.isNull());

   database->putInteger("PDAT_CELLDATA_VERSION", PDAT_CELLDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   array_database = database->putDatabase("d_data");
   d_data.putToDatabase(array_database);
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
