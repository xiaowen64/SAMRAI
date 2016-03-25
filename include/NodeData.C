//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated node centered patch data type
//

#ifndef included_pdat_NodeData_C
#define included_pdat_NodeData_C

#include "NodeData.h"

#include "Box.h"
#include "BoxList.h"
#include "NodeOverlap.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"

#define PDAT_NODEDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "NodeData.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for node data objects.  The constructor	*
* simply initializes data variables and sets up the array data.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
NodeData<DIM,TYPE>::NodeData(const hier::Box<DIM>& box,
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
   const hier::Box<DIM> node = NodeGeometry<DIM>::toNodeBox(this -> getGhostBox());
   d_data.initializeArray(node, depth, pool);

}

template<int DIM, class TYPE>
NodeData<DIM,TYPE>::~NodeData()
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
NodeData<DIM,TYPE>::NodeData(const NodeData<DIM,TYPE>& foo)
:  hier::PatchData<DIM>(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::operator=(const NodeData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Perform a fast copy between two node centered arrays where their	*
* index spaces overlap.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src)
{
   const NodeData<DIM,TYPE> *t_src =
      dynamic_cast<const NodeData<DIM,TYPE> *>(&src);
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
void NodeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst) const
{
   NodeData<DIM,TYPE> *t_dst =
      dynamic_cast<NodeData<DIM,TYPE> *>(&dst);
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
void NodeData<DIM,TYPE>::copy(const hier::PatchData<DIM>& src,
                                const hier::BoxOverlap<DIM>& overlap)
{
   const NodeData<DIM,TYPE> *t_src =
      dynamic_cast<const NodeData<DIM,TYPE> *>(&src);
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      d_data.copy(t_src->d_data,
                  t_overlap->getDestinationBoxList(),
                  t_overlap->getSourceOffset());
   }
}

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::copy2(hier::PatchData<DIM>& dst,
                                 const hier::BoxOverlap<DIM>& overlap) const
{
   NodeData<DIM,TYPE> *t_dst =
      dynamic_cast<NodeData<DIM,TYPE> *>(&dst);
   const NodeOverlap<DIM> *t_overlap = 
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
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
* Perform a fast copy from a node data object to this node data         *
* object at the specified depths, where their index spaces overlap.     *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::copyDepth(int dst_depth,
			           const NodeData<DIM,TYPE>& src,
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
bool NodeData<DIM,TYPE>::canEstimateStreamSizeFromBox() const
{
   return(ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox());
}

template<int DIM, class TYPE>
int NodeData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxOverlap<DIM>& overlap) const
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
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
void NodeData<DIM,TYPE>::packStream(tbox::AbstractStream& stream,
                                    const hier::BoxOverlap<DIM>& overlap) const
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   d_data.packStream(stream,
                     t_overlap->getDestinationBoxList(),
                     t_overlap->getSourceOffset());
}

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::unpackStream(tbox::AbstractStream& stream,
                                      const hier::BoxOverlap<DIM>& overlap)
{
   const NodeOverlap<DIM> *t_overlap =
      dynamic_cast<const NodeOverlap<DIM> *>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(t_overlap != NULL);
#endif
   d_data.unpackStream(stream,
                       t_overlap->getDestinationBoxList(),
                       t_overlap->getSourceOffset());
}

/*
*************************************************************************
*									*
* Calculate the amount of memory space needed to represent the data	*
* for a node centered grid.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t NodeData<DIM,TYPE>::getSizeOfData(const hier::Box<DIM>& box,
                                         int depth,
                                         const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   const hier::Box<DIM> ghost_box = hier::Box<DIM>::grow(box, ghosts);
   const hier::Box<DIM> node_box = NodeGeometry<DIM>::toNodeBox(ghost_box);
   return(ArrayData<DIM,TYPE>::getSizeOfData(node_box, depth));
}

/*
*************************************************************************
*                                                                       *
* Print node-centered data.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::print(const hier::Box<DIM>& box, 
                               std::ostream& os, 
                               int prec) const
{
   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      print(box, d, os, prec);
   }
}

template <int DIM, class TYPE>
void NodeData<DIM,TYPE>::print(const hier::Box<DIM>& box,
			       int depth,
			       std::ostream& os,
			       int prec) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   os.precision(prec);
   for (NodeIterator<DIM> i(box); i; i++) {
      os << "array" << i() << " = " << d_data(i(),depth) << std::endl << std::flush; 
   }
}

/*
*************************************************************************
*                                                                       *
* Checks to make sure that the class version and restart file		*
* version are equal.  If so, reads in d_depth and has d_data		*
* retrieve its own data from the database.				*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::getSpecializedFromDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("PDAT_NODEDATA_VERSION");
   if (ver != PDAT_NODEDATA_VERSION) {
      TBOX_ERROR("NodeData<DIM>::getSpecializedFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   tbox::Pointer<tbox::Database> array_database;
   array_database = database->getDatabase("d_data");
   (d_data).getFromDatabase(array_database);
}

/*
*************************************************************************
*                                                                       *
* Writes out the class version number and d_depth, Then has d_data	*
* write its own data to the database.					*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void NodeData<DIM,TYPE>::putSpecializedToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("PDAT_NODEDATA_VERSION", PDAT_NODEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   tbox::Pointer<tbox::Database> array_database;
   array_database = database->putDatabase("d_data");
   (d_data).putToDatabase(array_database);
}


}
}

#endif
