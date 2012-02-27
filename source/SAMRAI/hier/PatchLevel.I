/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A collection of patches at one level of the AMR hierarchy
 *
 ************************************************************************/

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getLevelNumber() const
{
   return d_level_number;
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getNextCoarserHierarchyLevelNumber() const
{
   return d_next_coarser_level_number;
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevel::inHierarchy() const
{
   return d_in_hierarchy;
}

SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<Patch>&
PatchLevel::getPatch(
   const GlobalId& gid) const
//   const BlockId& bid) const
{
   BoxId mbid(gid);
   PatchContainer::const_iterator it = d_patches.find(mbid);
   TBOX_ASSERT(it != d_patches.end());
   return it->second;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<Patch>
PatchLevel::getPatch(
   const BoxId& mapped_box_id) const
{
   const PatchContainer::const_iterator mi = d_patches.find(mapped_box_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   if (mi == d_patches.end()) {
      TBOX_ERROR("PatchLevel::getPatch(" << mapped_box_id
                                         << "): patch does not exist locally.");
   }
#endif
   return (*mi).second;
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getLocalNumberOfPatches() const
{
   return static_cast<int>(d_mapped_box_level->getLocalNumberOfBoxes());
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getGlobalNumberOfPatches() const
{
   return d_mapped_box_level->getGlobalNumberOfBoxes();
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getNumberOfPatches() const {
   return getGlobalNumberOfPatches();
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getLocalNumberOfCells() const
{
   return static_cast<int>(d_mapped_box_level->getLocalNumberOfCells());
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getGlobalNumberOfCells() const
{
   return d_mapped_box_level->getGlobalNumberOfCells();
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchDescriptor>
PatchLevel::getPatchDescriptor() const
{
   return d_descriptor;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchFactory>
PatchLevel::getPatchFactory() const
{
   return d_factory;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<GridGeometry>
PatchLevel::getGridGeometry() const
{
   return d_geometry;
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoxContainer>&
PatchLevel::getPhysicalDomainArray() const
{
   return d_physical_domain;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer&
PatchLevel::getPhysicalDomain(
   const BlockId& block_id) const
{
   return d_physical_domain[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
const IntVector&
PatchLevel::getRatioToLevelZero() const
{
   return d_ratio_to_level_zero;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
PatchLevel::getRatioToCoarserLevel() const
{
   return d_ratio_to_coarser_level;
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setRatioToCoarserLevel(
   const IntVector& ratio)
{
   d_ratio_to_coarser_level = ratio;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer&
PatchLevel::getBoxes() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedBoxLevel();
   }
   return d_boxes;
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::getBoxes(
   BoxContainer& boxes,
   const BlockId& block_id) const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedBoxLevel();
   }

   boxes.clear();
   const BoxContainer& global_mapped_boxes =
      d_mapped_box_level->getGlobalizedVersion().getGlobalBoxes();

   for (BoxContainerSingleBlockIterator gi(global_mapped_boxes, block_id);
        gi.isValid(); gi++) {
      boxes.pushBack(*gi);
   }
}

SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<BoxLevel>&
PatchLevel::getBoxLevel() const
{
   return d_mapped_box_level;
}

SAMRAI_INLINE_KEYWORD
const BoxLevel&
PatchLevel::getGlobalizedBoxLevel() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedBoxLevel();
   }
   return d_mapped_box_level->getGlobalizedVersion();
}

SAMRAI_INLINE_KEYWORD
const ProcessorMapping&
PatchLevel::getProcessorMapping() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedBoxLevel();
   }
   return d_mapping;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
PatchLevel::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
int
PatchLevel::getMappingForPatch(
   const BoxId& mapped_box_id) const
{
   // Note: p is required to be a local index.
   /*
    * This must be for backward compatability, because if p is a local
    * index, the mapping is always to d_mapped_box_level->getRank().
    * Here is the old code:
    *
    * return d_mapped_box_level->getBoxStrict(p)->getOwnerRank();
    */
   NULL_USE(mapped_box_id);
   return d_mapped_box_level->getMPI().getRank();
}

SAMRAI_INLINE_KEYWORD
const Box&
PatchLevel::getBoxForPatch(
   const BoxId& mapped_box_id) const
{
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mapped_box_level->getMPI().getRank());
   return getPatch(mapped_box_id)->getBox();
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevel::patchTouchesRegularBoundary(
   const BoxId& mapped_box_id) const
{
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mapped_box_level->getMPI().getRank());
   bool rval = getPatch(mapped_box_id)->getPatchGeometry()->getTouchesRegularBoundary();
   return rval;
}

/*
 * ************************************************************************
 *
 * Allocate or deallocate data for single components or collections of
 * component on all patches on a patch level.
 *
 * ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::allocatePatchData(
   const int id,
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->allocatePatchData(id, timestamp);
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::allocatePatchData(
   const ComponentSelector& components,
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->allocatePatchData(components, timestamp);
   }
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevel::checkAllocated(
   const int id) const
{
   bool allocated = true;
   for (PatchContainer::const_iterator
        mi = d_patches.begin(); mi != d_patches.end(); ++mi) {
      allocated &= (*mi).second->checkAllocated(id);
   }
   return allocated;
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::deallocatePatchData(
   const int id)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->deallocatePatchData(id);
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::deallocatePatchData(
   const ComponentSelector& components)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->deallocatePatchData(components);
   }
}

/*
 * ************************************************************************
 *
 * Set the simulation time for all patches in the patch level.
 *
 * ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setTime(
   const double timestamp,
   const int id)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp, id);
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setTime(
   const double timestamp,
   const ComponentSelector& components)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp, components);
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setTime(
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp);
   }
}

/*
 * ************************************************************************
 *
 * Set level numbers relating this level to "level", a level in
 * a hierarchy
 *
 * ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setLevelNumber(
   const int level)
{
   d_level_number = level;

   for (PatchLevel::Iterator p(this); p; p++) {
      p->setPatchLevelNumber(d_level_number);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setNextCoarserHierarchyLevelNumber(
   const int level)
{
   d_next_coarser_level_number = level;
}

/*
 * ************************************************************************
 *
 * Set whether this level resides in a hierarchy.
 *
 * ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setLevelInHierarchy(
   bool in_hierarchy)
{
   d_in_hierarchy = in_hierarchy;

   for (PatchLevel::Iterator p(this); p; p++) {
      p->setPatchInHierarchy(d_in_hierarchy);
   }
}

/*
 * ************************************************************************
 *
 * Call the geometry routine to create and set boundary boxes, if they
 * have not already been created.
 *
 * ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
PatchLevel::setBoundaryBoxes()
{
   if (!d_boundary_boxes_created) {
      d_geometry->setBoundaryBoxes(*this);
      d_boundary_boxes_created = true;
   }
}

// Iterator implementation

/*
 *************************************************************************
 * Assignment.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator&
PatchLevel::Iterator::operator = (
   const PatchLevel::Iterator& rhs)
{
   d_iterator = rhs.d_iterator;
   d_patches = rhs.d_patches;
   return *this;
}

/*
 *************************************************************************
 * Standard iterator's dereference operator.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<Patch>&
PatchLevel::Iterator::operator * () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Alternative iterator's dereference operator.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<Patch>&
PatchLevel::Iterator::operator () () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Delegate boost::shared_ptr operations to the Patch pointer.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<Patch>&
PatchLevel::Iterator::operator -> () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Equality comparison.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
PatchLevel::Iterator::operator == (
   const Iterator& rhs) const
{
   return d_iterator == rhs.d_iterator;
}

/*
 *************************************************************************
 * Inequality comparison.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
PatchLevel::Iterator::operator != (
   const Iterator& rhs) const
{
   return d_iterator != rhs.d_iterator;
}

/*
 *************************************************************************
 * Standard iterator's pre-increment.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const PatchLevel::Iterator&
PatchLevel::Iterator::operator ++ ()
{
   ++d_iterator;
   return *this;
}

/*
 *************************************************************************
 * Standard iterator's post-increment.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator
PatchLevel::Iterator::operator ++ (
   int)
{
   Iterator tmp_iterator = *this;
   ++d_iterator;
   return tmp_iterator;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void
PatchLevel::Iterator::initialize(
   const boost::shared_ptr<PatchLevel>& patch_level)
{
   d_iterator = patch_level->d_patches.begin();
   d_patches = &patch_level->d_patches;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void
PatchLevel::Iterator::initialize(
   const PatchLevel& patch_level)
{
   d_iterator = patch_level.d_patches.begin();
   d_patches = &patch_level.d_patches;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void
PatchLevel::Iterator::initialize(
   const PatchLevel* patch_level)
{
   d_iterator = patch_level->d_patches.begin();
   d_patches = &patch_level->d_patches;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::operator bool ()
{
   return d_iterator != d_patches->end();
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
