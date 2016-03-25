//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/multiblock/MultiblockPatchLevel.h $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: class to manage multiblock levels
//

#ifndef included_hier_MultiblockPatchLevel
#define included_hier_MultiblockPatchLevel

#include "SAMRAI_config.h"
#include "BasePatchLevel.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/DescribedClass.h"


namespace SAMRAI {
    namespace hier {

/*!
 * @brief Class MultiblockPatchLevel<DIM> contains an array of
 * hier::PatchLevel<DIM> that contains all of the patch levels that have the
 * same level of refinement in a multiblock domain.
 *
 * @see hier::PatchLevel
 * @see hier::MultiblockPatchHierarchy
 */

template <int DIM>
class MultiblockPatchLevel
 : public hier::BasePatchLevel<DIM>
{

public:

   /*!
    * @brief Constructor takes an array of pointers to patch levels.
    *
    * @param levels Array of pointers to hier::PatchLevel<DIM>.  The array
    *               indices correspond to block numbers.  Pointers in the
    *               array may be null.  A null pointer indicates that the
    *               MultiblockPatchLevel does not represent any space in
    *               the block associated with its array index. 
    */
   MultiblockPatchLevel(
      tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >& levels);

   /*!
    * Destructor is uninteresting
    */
   ~MultiblockPatchLevel<DIM>();

   /*!
    * @brief Return a pointer to the hier::PatchLevel<DIM> associated with the
    * id number.
    *
    * @param id corresponds to the indexing of the array passed into the
    *           constructor.
    */
   tbox::Pointer< hier::PatchLevel<DIM> >
   getPatchLevelForBlock(const int id) const;

   /*!
    * @brief Allocate the specified component on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    *
    * @param id        A patch data id
    * @param timestamp Simulation time
    * @param pool      A memory arena
    */
   void allocatePatchData(
      const int id,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL);

   /*!
    * @brief Allocate the specified components on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    *
    * @param components A component selector defining a set of patch data id's
    * @param timestamp  Simulation time
    * @param pool       A memory arena
    */
   void allocatePatchData(
      const hier::ComponentSelector& components,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL);

   /*!
    * @brief Deallocate the specified component on all patches
    *
    * @param id Patch data id of data to be deallocated
    */ 
   void deallocatePatchData(const int id);

   /*!
    * @brief Deallocate the specified components on all patches
    *
    * @param components A component selector defining a set of patch data id's
    *                   which indicate the data to be deallocated
    */ 
   void deallocatePatchData(const hier::ComponentSelector& components);

   /*!
    * @brief Set the simulation time for the specified patch component.
    *
    * @param timestamp Simulation time
    * @param id        A patch data id
    */
   void setTime(const double timestamp, const int id);

   /*!
    * @brief Set the simulation time for the specified patch components.
    *
    * @param timestamp  Simulation time
    * @param components A component selector defining a set of patch data id's
    */
   void setTime(
      const double timestamp,
      const hier::ComponentSelector& components);

   /*!
    * @brief Set the simulation time for all allocated patch components.
    *
    * @param timestamp Simulation time
    */
   void setTime(const double timestamp);

   /*!
    * @brief Get the number of blocks in the multiblock domain
    */
   int getNumberOfBlocks() const;

   /*!
    * @brief Get the level number of this level
    */
   int getLevelNumber() const;

   /*!
    * @brief Get the ratio to level zero of this level
    */
   const hier::IntVector<DIM>& getRatio() const;

private:

   /*
    * Array of hierarchies that cover the multiblock domain
    */
   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > > d_levels;

   int d_number_blocks;
};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockPatchLevel.C"
#endif

#endif
