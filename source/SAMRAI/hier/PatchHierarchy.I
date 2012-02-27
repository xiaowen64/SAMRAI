/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An AMR hierarchy of patch levels
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchLevel>
PatchHierarchy::getPatchLevel(
   const int l) const
{
   TBOX_ASSERT((l >= 0) && (l < d_number_levels));

   return d_patch_levels[l];
}

SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<BoxLevel>&
PatchHierarchy::getBoxLevel(
   int level) const
{
   return d_patch_levels[level]->getBoxLevel();
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchDescriptor>
PatchHierarchy::getPatchDescriptor() const
{
   return d_patch_descriptor;
}

SAMRAI_INLINE_KEYWORD
const BoxLevel&
PatchHierarchy::getDomainBoxLevel() const
{
   return d_domain_mapped_box_level;
}

SAMRAI_INLINE_KEYWORD
const tbox::SAMRAI_MPI&
PatchHierarchy::getMPI() const
{
   return d_domain_mapped_box_level.getMPI();
}

SAMRAI_INLINE_KEYWORD
bool
PatchHierarchy::levelExists(
   const int l) const
{
   return (l < d_number_levels) && d_patch_levels[l];
}

SAMRAI_INLINE_KEYWORD
bool
PatchHierarchy::finerLevelExists(
   const int l) const
{
   return (l + 1 < d_number_levels) && d_patch_levels[l + 1];
}

SAMRAI_INLINE_KEYWORD
bool
PatchHierarchy::levelCanBeRefined(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   return level_number < getMaxNumberOfLevels() - 1;
}

SAMRAI_INLINE_KEYWORD
int
PatchHierarchy::getNumberOfLevels() const
{
   return d_number_levels;
}

SAMRAI_INLINE_KEYWORD
int
PatchHierarchy::getMaxNumberOfLevels() const
{
   return d_max_levels;
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setMaxNumberOfLevels(
   int max_levels)
{
   d_max_levels = max_levels;

   if (d_max_levels != int(d_ratio_to_coarser.size())) {
      d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
      d_smallest_patch_size.resize(d_max_levels, d_smallest_patch_size.back());
      d_largest_patch_size.resize(d_max_levels, d_largest_patch_size.back());
      d_proper_nesting_buffer.resize(d_max_levels - 1,
         d_proper_nesting_buffer.empty() ? 1 : d_proper_nesting_buffer.back());
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setRatioToCoarserLevel(
   const IntVector& ratio,
   int level)
{
   TBOX_ASSERT(level > 0 && level < d_max_levels);
   d_ratio_to_coarser[level] = ratio;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
PatchHierarchy::getRatioToCoarserLevel(
   int level) const
{
   TBOX_ASSERT(level < d_max_levels);
   return d_ratio_to_coarser[level];
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setSmallestPatchSize(
   const IntVector& size,
   int level)
{
   TBOX_ASSERT(level >= 0 && level < d_max_levels);
   d_smallest_patch_size[level] = size;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
PatchHierarchy::getSmallestPatchSize(
   int level) const
{
   TBOX_ASSERT(level >= 0 && level < d_max_levels);
   return d_smallest_patch_size[level];
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setLargestPatchSize(
   const IntVector& size,
   int level)
{
   TBOX_ASSERT(level >= 0 && level < d_max_levels);
   d_largest_patch_size[level] = size;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
PatchHierarchy::getLargestPatchSize(
   int level) const
{
   TBOX_ASSERT(level >= 0 && level < d_max_levels);
   return d_largest_patch_size[level];
}

SAMRAI_INLINE_KEYWORD
int
PatchHierarchy::getProperNestingBuffer(
   int ln) const
{
   TBOX_ASSERT(ln >= 0 && ln < d_max_levels);
   return (ln < d_max_levels - 1) ? d_proper_nesting_buffer[ln] : -1;
}

SAMRAI_INLINE_KEYWORD
bool
PatchHierarchy::allowPatchesSmallerThanGhostWidth() const
{
   return d_allow_patches_smaller_than_ghostwidth;
}

SAMRAI_INLINE_KEYWORD
bool
PatchHierarchy::allowPatchesSmallerThanMinimumSize() const
{
   return d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps;
}

SAMRAI_INLINE_KEYWORD
int
PatchHierarchy::getFinestLevelNumber() const
{
   return d_number_levels - 1;
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setPatchFactory(
   const boost::shared_ptr<PatchFactory>& factory)
{
   d_patch_factory = factory;
}

SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::setPatchLevelFactory(
   const boost::shared_ptr<PatchLevelFactory>& factory)
{
   d_patch_level_factory = factory;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<GridGeometry>
PatchHierarchy::getGridGeometry() const
{
   return d_grid_geometry;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
PatchHierarchy::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
const std::string&
PatchHierarchy::getObjectName() const
{
   return d_object_name;
}


SAMRAI_INLINE_KEYWORD
void
PatchHierarchy::initializeCallback()
{
   /*
    * No-op.  This class doesn't
    */
}

}
}
