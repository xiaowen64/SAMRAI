/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2015 Lawrence Livermore National Security, LLC
 * Description:   Manages PatchData objects created at coarse-fine stencils
 *
 ************************************************************************/

#ifndef included_xfer_CompositeBoundarySchedule
#define included_xfer_CompositeBoundarySchedule

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class CompositeBoundarySchedule manages the communication and access
 * of patch data objects that hold data at a composite boundary.
 *
 * The purpose of this class is to allow an application code that is operating
 * on a particular Patch from a PatchHierarchy to have a local view into data
 * from the next finer level that exists within a certain stencil width of the
 * Patch's coarse-fine boundaries (if it has any).
 *
 * The recommended usage is to create a CompositeBoundaryAlgorithm to hold
 * the stencil width and the patch data ids, and use
 * CompositeBoundaryAlgorithm::createSchedule to create instances of this
 * class.  
 */

class CompositeBoundarySchedule
{
public:

   /*!
    * @brief Constructor sets up basic state
    *
    * The constructor associates the CompositeBoundarySchedule with a hierarchy
    * and a stencil width.  It does not manage any data until further
    * setup calls are made.
    *
    * The stencil width is interpreted in terms of the resolution of finer
    * levels.  So if the value is 4, that means that code operating on level n
    * will be able to see a coarse-fine stencil of width 4 in the mesh 
    * resolution of level n + 1;
    *
    * This object will become invalid if either the coarse level or the
    * next finer level is changed by a regrid of the hierarchy.
    *
    * @param[in] hierarchy
    * @param[in] coarse_level_num  Level number of the coarse level at
    *                              coarse-fine boundary
    * @param[in] stencil_width  Width in terms of fine level zones
    * @param[in] data_ids   Patch data ids for data to be communicated 
    *
    * @pre hierarchy != 0
    * @pre stencil_width >= 1
    */
   CompositeBoundarySchedule(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int coarse_level_num,
      int stencil_width,
      const std::set<int>& data_ids);

   /*!
    * @brief Destructor
    */
   ~CompositeBoundarySchedule();

   /*!
    * @brief Fill data into the stencil that exists at the coarse-fine
    * boundaries of the given level.
    *
    * Data willl be filled from the corresponding fine level from the
    * onto the stencil, where it is then accessible through the method
    * getBoundaryPatchData().
    *
    * This may be called multiple times, to account for changing state of
    * data on the hierarchy.
    */
   void fillData();

   /*!
    * @brief Get PatchData from the composite boundary stencil for a given
    *        data id.
    *
    * Get a vector that holds the stencil PatchData for the parts of the
    * coarse-fine boundary touching the Patch.  The Patch must be a local
    * Patch on the hierarchy, and fillData() must already have been called on
    * this object.
    *
    * If the Patch does not touch any coarse-fine boundaries, the returned
    * vector will be empty.
    *
    * @param patch    a local Patch from the hierarchy
    * @param data_id  patch data id for data stored in the stencil
    */  
   const std::vector<boost::shared_ptr<hier::PatchData> >&
   getBoundaryPatchData(
      const hier::Patch& patch,
      int data_id) const;

private:
   CompositeBoundarySchedule(
      const CompositeBoundarySchedule&);                  // not implemented
   CompositeBoundarySchedule&
   operator = (
      const CompositeBoundarySchedule&);                  // not implemented

   /*!
    * @brief Add a patch data id to identify data to be managed on the stencils.
    *
    * This may be called multiple times to add as many data ids as desired.
    * All desired data ids must be added prior to calling
    * createStencilForLevel().
    *
    * @param[in] data_id
    *
    * @pre data_id >= 0 
    */
   void addDataId(int data_id);

   /*!
    * @brief Create a stencil that exists at the coarse-fine boundaries
    * of the given level.
    *
    * The given level number represents a level from the hierarchy.  The
    * stencil will hold data from the next finer level of resolution at the
    * coarse fine boundaries.
    *
    * This method should be called to recreate stencils after the level
    * represented by level_num and/or the next finer level have been regridded.
    *
    * @param level_num   level number for a level from hierarchy.
    */
   void createStencilForLevel();

   /*!
    * @brief Set up communication schedule to fill data at composite boundary
    *
    * This sets up the communications that fill the patch data around the
    * coarse-fine boundaries of the given level.
    *
    * @param[in]   level_num
    *
    */
   void setUpSchedule();


   /*!
    * @brief Hierarchy where the stencils exist
    */
   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Levels representing the stencils
    */
   boost::shared_ptr<hier::PatchLevel> d_stencil_level;

   /*!
    * @brief Map from BoxId of a hierarchy patch to the BoxIds of its
    * coarse-fine stencil boxes.
    */
   std::map<hier::BoxId, std::set<hier::BoxId> > d_stencil_map;

   /*!
    * @brief Width of the stencils.  The width represents a width in
    * terms of fine resolution at coarse-fine boundaries.
    */
   int d_stencil_width;

   /*!
    * @brief Patch data ids for data to be held on the stencil.
    */
   std::set<int> d_data_ids;

   /*!
    * @brief Algorithm for communication of data to the stencil.
    */
   RefineAlgorithm d_refine_algorithm;

   /*!
    * @brief Schedules for communication of data from hierarchy to stencils.
    */
   boost::shared_ptr<RefineSchedule> d_refine_schedule;

   /*
    * typedefs to simplify nested container declaration.
    */
   typedef std::vector<boost::shared_ptr<hier::PatchData> > VectorPatchData;
   typedef std::map<hier::BoxId, VectorPatchData> MapBoxIdPatchData;

   /*!
    * @brief Container of PatchData on the stencil.
    *
    * The inner nested container is a vector of PatchData which consists of 
    * one PatchData object for each box of the coarse-fine stencil for a single
    * Patch of the hierarchy.
    *
    * The nesting of the containers is organized as:
    *
    * d_data_map[level_number][data_id][box_id][stencil_box_index]
    *
    * level_number is the level number of the coarse level at a coarse-fine
    * boundary where the stencil has been created and filled.  data_id
    * is the patch data id of the data that is to be accessed.  box_id is
    * the id of a Patch on the hierarchy.  stencil_box_index is an index into
    * the inner nested vector of PatchData.
    */
   std::map< int, MapBoxIdPatchData > d_data_map;

   /*!
    * @brief Vector of flags telling if stencil has been created for each
    * level of the hierarchy.
    */
   bool d_stencil_created;

   /*!
    * @brief Vector of flags telling if stencil has been filled for each
    * level of the hierarchy.
    */
   bool d_stencil_filled; 

   int d_coarse_level_num;
};

}
}

#endif
