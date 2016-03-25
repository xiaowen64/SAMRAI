//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockGridGeometry.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description:	GridGeometry for Multiblock.
//
 
#ifndef included_hier_MultiblockGridGeometry
#define included_hier_MultiblockGridGeometry

#include "SAMRAI_config.h"

#include "GridGeometry.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace hier {

/*!
 * @brief Class MultiblockGridGeometry<DIM>
 */

template<int DIM>
class MultiblockGridGeometry
{
public:
   /*!
    * @brief The constructor takes an array of GridGeometry pointers
    */
   MultiblockGridGeometry(
      tbox::Array< tbox::Pointer< hier::GridGeometry<DIM> > >& block_geoms);
 
   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockGridGeometry<DIM>();

   tbox::Pointer< hier::GridGeometry<DIM> >& getBlockGeometry(const int block)
   {
      return d_block_geometry[block];
   }

   tbox::Array< tbox::Pointer< hier::GridGeometry<DIM> > >&
   getBlockGeometryArray()
   {
      return d_block_geometry;
   }

private:

   tbox::Array< tbox::Pointer< hier::GridGeometry<DIM> > > d_block_geometry;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockGridGeometry.C"
#endif

#endif
