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
int PatchLevel::getLevelNumber() const
{
   return d_level_number;
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getNextCoarserHierarchyLevelNumber() const
{
   return d_next_coarser_level_number;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevel::inHierarchy() const
{
   return d_in_hierarchy;
}

SAMRAI_INLINE_KEYWORD
const tbox::Pointer<Patch>& PatchLevel::getPatch(
   const GlobalId& gid,
   const BlockId& bid) const
{
   BoxId mbid(gid, bid);
   PatchContainer::const_iterator it = d_patches.find(mbid);
   TBOX_ASSERT(it != d_patches.end());
   return it->second;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<Patch> PatchLevel::getPatch(
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
int PatchLevel::getLocalNumberOfPatches() const
{
   return static_cast<int>(d_mapped_box_level->getLocalNumberOfBoxes());
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getGlobalNumberOfPatches() const
{
   return d_mapped_box_level->getGlobalNumberOfBoxes();
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getNumberOfPatches() const {
   return getGlobalNumberOfPatches();
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getLocalNumberOfCells() const
{
   return static_cast<int>(d_mapped_box_level->getLocalNumberOfCells());
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getGlobalNumberOfCells() const
{
   return d_mapped_box_level->getGlobalNumberOfCells();
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchDescriptor> PatchLevel::getPatchDescriptor() const
{
   return d_descriptor;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchFactory> PatchLevel::getPatchFactory() const
{
   return d_factory;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<GridGeometry> PatchLevel::getGridGeometry() const
{
   return d_geometry;
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoxList>& PatchLevel::getPhysicalDomainArray() const
{
   return d_physical_domain;
}

SAMRAI_INLINE_KEYWORD
const BoxList& PatchLevel::getPhysicalDomain(const BlockId& block_id) const
{
   return d_physical_domain[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
const IntVector& PatchLevel::getRatioToLevelZero() const
{
   return d_ratio_to_level_zero;
}

SAMRAI_INLINE_KEYWORD
const IntVector& PatchLevel::getRatioToCoarserLevel() const
{
   return d_ratio_to_coarser_level;
}

SAMRAI_INLINE_KEYWORD
void PatchLevel::setRatioToCoarserLevel(
   const IntVector& ratio)
{
   d_ratio_to_coarser_level = ratio;
}

SAMRAI_INLINE_KEYWORD
const BoxList& PatchLevel::getBoxes() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedMappedBoxLevel();
   }
   return d_boxes;
}

SAMRAI_INLINE_KEYWORD
void PatchLevel::getBoxes(BoxList& boxes, const BlockId& block_id) const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedMappedBoxLevel();
   }

   boxes.clearItems();
   const BoxSet& global_mapped_boxes =
      d_mapped_box_level->getGlobalizedVersion().getGlobalMappedBoxes();

   for (BoxSetSingleBlockIterator gi(global_mapped_boxes, block_id);
        gi.isValid(); gi++) {
      boxes.appendItem(*gi);
   }
}

SAMRAI_INLINE_KEYWORD
const tbox::ConstPointer<MappedBoxLevel>& PatchLevel::getMappedBoxLevel() const
{
   return d_mapped_box_level;
}

SAMRAI_INLINE_KEYWORD
const MappedBoxLevel&
PatchLevel::getGlobalizedMappedBoxLevel() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedMappedBoxLevel();
   }
   return d_mapped_box_level->getGlobalizedVersion();
}

#if 0
// Removed due to incompatibility with DLBG.
SAMRAI_INLINE_KEYWORD
const tbox::Array<tbox::List<IntVector> >& PatchLevel::getShiftsForLevel()
const
{
   return d_shifts;
}
#endif

SAMRAI_INLINE_KEYWORD
const ProcessorMapping& PatchLevel::getProcessorMapping() const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedMappedBoxLevel();
   }
   return d_mapping;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& PatchLevel::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
int PatchLevel::getMappingForPatch(
   const BoxId& mapped_box_id) const
{
   // Note: p is required to be a local index.
   /*
    * This must be for backward compatability, because if p is a local
    * index, the mapping is always to d_mapped_box_level->getRank().
    * Here is the old code:
    *
    * return d_mapped_box_level->getMappedBoxStrict(p)->getOwnerRank();
    */
   NULL_USE(mapped_box_id);
   return d_mapped_box_level->getRank();
}

SAMRAI_INLINE_KEYWORD
const Box& PatchLevel::getBoxForPatch(
   const BoxId& mapped_box_id) const
{
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mapped_box_level->getRank());
   return getPatch(mapped_box_id)->getBox();
}

SAMRAI_INLINE_KEYWORD
bool PatchLevel::patchTouchesRegularBoundary(
   const BoxId& mapped_box_id) const
{
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mapped_box_level->getRank());
   bool rval = getPatch(mapped_box_id)->getPatchGeometry()->getTouchesRegularBoundary();
   return rval;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevel::patchTouchesPeriodicBoundary(
   const BoxId& mapped_box_id) const
{
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mapped_box_level->getRank());
   bool rval = getPatch(mapped_box_id)->getPatchGeometry()->getTouchesPeriodicBoundary();
   return rval;
}

// Iterator implementation

/*
 *************************************************************************
 * Default constructor.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator():
   d_iterator(),
   d_patches(NULL /* Unused since not backward compatibility not needed */)
{
}

/*
 *************************************************************************
 * Copy constructor.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator(
   const PatchLevel::Iterator& r):
   d_iterator(r.d_iterator),
   d_patches(NULL /* Unused since not backward compatibility not needed */)
{
}

/*
 *************************************************************************
 * Construct from raw iterator.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator(
   const PatchLevel::PatchContainer::const_iterator& r):
   d_iterator(r),
   d_patches(NULL /* Unused since not backward compatibility not needed */)
{
}

/*
 *************************************************************************
 * Assignment.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator& PatchLevel::Iterator::operator = (
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
const tbox::Pointer<Patch>& PatchLevel::Iterator::operator * () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Alternative iterator's dereference operator.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Pointer<Patch>& PatchLevel::Iterator::operator () () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Delegate Pointer operations to the Patch pointer.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Pointer<Patch>& PatchLevel::Iterator::operator -> () const
{
   return d_iterator->second;
}

/*
 *************************************************************************
 * Equality comparison.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PatchLevel::Iterator::operator == (
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
bool PatchLevel::Iterator::operator != (
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
const PatchLevel::Iterator& PatchLevel::Iterator::operator ++ ()
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
PatchLevel::Iterator PatchLevel::Iterator::operator ++ (
   int)
{
   Iterator tmp_iterator = d_iterator;
   ++d_iterator;
   return tmp_iterator;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator(
   const PatchLevel& patch_level):
   d_iterator(patch_level.d_patches.begin()),
   d_patches(&patch_level.d_patches)
{
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator(
   const tbox::Pointer<PatchLevel>& patch_level):
   d_iterator(patch_level->d_patches.begin()),
   d_patches(&patch_level->d_patches)
{
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
PatchLevel::Iterator::Iterator(
   const PatchLevel* patch_level):
   d_iterator(patch_level->d_patches.begin()),
   d_patches(&patch_level->d_patches)
{
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void PatchLevel::Iterator::initialize(
   const tbox::Pointer<PatchLevel>& patch_level)
{
   d_iterator = patch_level->d_patches.begin();
   d_patches = &patch_level->d_patches;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void PatchLevel::Iterator::initialize(
   const PatchLevel& patch_level)
{
   d_iterator = patch_level.d_patches.begin();
   d_patches = &patch_level.d_patches;
}

// Support for backward compatible interface by new PatchLevel::Iterator.
SAMRAI_INLINE_KEYWORD
void PatchLevel::Iterator::initialize(
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
