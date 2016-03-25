//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/multiblock/MultiblockPatchLevel.C $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Base class for geometry management on patches
//

#ifndef included_hier_MultiblockPatchLevel_C
#define included_hier_MultiblockPatchLevel_C

#include "MultiblockPatchLevel.h"


namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for multiblcok level.  The constructor     *
* simply initializes the d_levels data member.                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> MultiblockPatchLevel<DIM>::MultiblockPatchLevel(
   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >& levels)
{

   d_levels = levels;
   d_number_blocks = levels.getSize();
}

template<int DIM> MultiblockPatchLevel<DIM>::~MultiblockPatchLevel()
{
}

/*
*************************************************************************
*                                                                       *
* Returns a single patch level from the block identified by the integer *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::PatchLevel<DIM> >
MultiblockPatchLevel<DIM>::getPatchLevelForBlock(const int id) const
{
   return d_levels[id];
}

/*
*************************************************************************
*                                                                       *
* Returns the number of blocks.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> int MultiblockPatchLevel<DIM>::getNumberOfBlocks() const
{
   return (d_number_blocks);
}

/*
*************************************************************************
*                                                                       *
* Returns level number.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> int MultiblockPatchLevel<DIM>::getLevelNumber() const
{
   int level_num = 0;
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         level_num = d_levels[i]->getLevelNumber();
         break;
      }
   }
   return (level_num);
}

/*
*************************************************************************
*                                                                       *
* Allocate patch data for the given id                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MultiblockPatchLevel<DIM>::allocatePatchData(
   const int id,
   const double timestamp,
   tbox::Pointer<tbox::Arena> pool)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->allocatePatchData(id, timestamp, pool);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Allocate patch data for the given components                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MultiblockPatchLevel<DIM>::allocatePatchData(
   const hier::ComponentSelector& components,
   const double timestamp,
   tbox::Pointer<tbox::Arena> pool)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->allocatePatchData(components, timestamp, pool);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Deallocate patch data for the given id                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockPatchLevel<DIM>::deallocatePatchData(const int id)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->deallocatePatchData(id);
      }
   }
}
 
/*
*************************************************************************
*                                                                       *
* Deallocate patch data for the given components                        *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockPatchLevel<DIM>::deallocatePatchData(
   const hier::ComponentSelector& components)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->deallocatePatchData(components);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set simulation time                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockPatchLevel<DIM>::setTime(const double timestamp, const int id)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->setTime(timestamp, id);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set simulation time                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockPatchLevel<DIM>::setTime(const double timestamp,
                                    const hier::ComponentSelector& components)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->setTime(timestamp, components);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set simulation time                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockPatchLevel<DIM>::setTime(const double timestamp)
{
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         d_levels[i]->setTime(timestamp);
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Get ratio to level zero                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
const hier::IntVector<DIM>& MultiblockPatchLevel<DIM>::getRatio() const
{
   int block_num = 0;
   for (int i=0; i < d_levels.getSize(); i++) {
      if (!d_levels[i].isNull()) {
         block_num = i;
         break;
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_levels[block_num].isNull());
#endif

   return(d_levels[block_num]->getRatio());
}

}
}
#endif
