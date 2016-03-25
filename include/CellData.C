//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/cell/CellData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated cell centered patch data type
//

#ifndef included_pdat_CellData_C
#define included_pdat_CellData_C

#include "CellData.h"

#include "Box.h"
#include "BoxList.h"
#include "CellGeometry.h"
#include "CellOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"

#define PDAT_CELLDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "CellData.I"
#endif


namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Calculate the amount of memory space needed to represent the data	*
* for a cell centered grid.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t CellData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box,
   int depth,
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   const hier::Box<DIM> ghost_box = hier::Box<DIM>::grow(box, ghosts);
   return(ArrayData<DIM,TYPE>::getSizeOfData(ghost_box, depth));
}

/*
*************************************************************************
*									*
* Constructor and destructor for cell data objects.  The constructor	*
* simply initializes data variables and sets up the array data.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
CellData<DIM,TYPE>::CellData(const hier::Box<DIM>& box,
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
   d_data.initializeArray(this -> getGhostBox(), depth, pool);
}

template<int DIM, class TYPE>
CellData<DIM,TYPE>::~CellData()
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
CellData<DIM,TYPE>::CellData(const CellData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::operator=(const CellData<DIM,TYPE>& foo)
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

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const CellData<DIM,TYPE> *t_src =
      dynamic_cast<const CellData<DIM,TYPE> *>(&src);
   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      const hier::Box<DIM> box = d_data.getBox() * t_src->d_data.getBox();
      if (!box.empty()) {
         d_data.copy(t_src->d_data, box);
      }
   }
}

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   CellData<DIM,TYPE> *t_dst = dynamic_cast<CellData<DIM,TYPE> *>(&dst);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
#endif
   const hier::Box<DIM> box = d_data.getBox() * t_dst->d_data.getBox();
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

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::copy(
   const hier::PatchData<DIM>& src, const hier::BoxOverlap<DIM>& overlap)
{
   const CellData<DIM,TYPE> *t_src = 
      dynamic_cast<const CellData<DIM,TYPE> *>(&src);

   const CellOverlap<DIM> *t_overlap = 
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);

   if ( (t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      d_data.copy(t_src->d_data,
                  t_overlap->getDestinationBoxList(),
                  t_overlap->getSourceOffset());
   }
}

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::copy2(
   hier::PatchData<DIM>& dst, const hier::BoxOverlap<DIM>& overlap) const
{
   CellData<DIM,TYPE> *t_dst = dynamic_cast<CellData<DIM,TYPE> *>(&dst);
   const CellOverlap<DIM> *t_overlap = 
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);
#endif
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

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::copyDepth(int dst_depth,
				     const CellData<DIM,TYPE>& src,
				     int src_depth)
{
  const hier::Box<DIM> box = d_data.getBox() * src.d_data.getBox();
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

template<int DIM, class TYPE>
bool CellData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int CellData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const CellOverlap<DIM> *t_overlap = 
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   return( d_data.getDataStreamSize(t_overlap->getDestinationBoxList(),
                                    t_overlap->getSourceOffset()) );
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
void CellData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream, const hier::BoxOverlap<DIM>& overlap) const
{
   const CellOverlap<DIM> *t_overlap = 
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   d_data.packStream( stream, t_overlap->getDestinationBoxList(), 
                              t_overlap->getSourceOffset());
}

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream, const hier::BoxOverlap<DIM>& overlap)
{
   const CellOverlap<DIM> *t_overlap = 
      dynamic_cast<const CellOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   d_data.unpackStream( stream, t_overlap->getDestinationBoxList(),
                                t_overlap->getSourceOffset() );
}

/*
*************************************************************************
*									*
* Print cell centered data.  Note:  makes call to specialized print     *
* routine in CellDataSpecialized.C       	 			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::print(
   const hier::Box<DIM>& box, 
   std::ostream& os, 
   int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      print(box, d, os, prec);
   }
}

template<int DIM, class TYPE> void CellData<DIM,TYPE>::print(
   const hier::Box<DIM>& box,
   int depth,
   std::ostream& os, 
   int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   os.precision(prec);
   for (CellIterator<DIM> i(box); i; i++) {
      os << "array" << i() << " = " << d_data(i(),depth) << std::endl << std::flush;
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

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

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

template<int DIM, class TYPE>
void CellData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_CELLDATA_VERSION", PDAT_CELLDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   array_database = database->putDatabase("d_data");
   d_data.putToDatabase(array_database);
}


}
}

#endif
