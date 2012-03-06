/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Spatial Key used for generating space-filling curves.
 *
 ************************************************************************/

namespace SAMRAI {
namespace mesh {

/*
 ****************************************************************************
 *
 * Copies the spatial key.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
SpatialKey&
SpatialKey::operator = (
   const SpatialKey& spatial_key)
{
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      d_key[i] = spatial_key.d_key[i];
   }
   return *this;
}

/*
 ****************************************************************************
 *
 * Compare two spatial keys for equality.  All values in d_key must be
 * equal for this to return true.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
SpatialKey::operator == (
   const SpatialKey& spatial_key) const
{
   bool are_equal = true;
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      if (d_key[i] != spatial_key.d_key[i]) {
         are_equal = false;
         break;
      }
   }
   return are_equal;
}

/*
 ****************************************************************************
 *
 * Compare two spatial keys for inequality.  Returns true if the keys are
 * not equal.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
SpatialKey::operator != (
   const SpatialKey& spatial_key) const
{
   return !((*this) == spatial_key);
}

/*
 ****************************************************************************
 *
 * Zeroes key by setting all values in d_key to 0.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
SpatialKey::setToZero()
{
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      d_key[i] = 0;
   }
}

/*
 ****************************************************************************
 *
 * Less than or equal operator for spatial keys.  Returns true if
 * this key is less or equal to the argument key.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
SpatialKey::operator <= (
   const SpatialKey& spatial_key) const
{
   return ((*this) < spatial_key) || ((*this) == spatial_key);
}

/*
 ****************************************************************************
 *
 * Greater than operator for spatial keys. Returns true if this key is
 * is greater than the argument key.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
SpatialKey::operator > (
   const SpatialKey& spatial_key) const
{
   return !((*this) < spatial_key) && ((*this) != spatial_key);
}

/*
 ****************************************************************************
 *
 * Greater than or equal operator for spatial keys.  Returns true if this
 * key is greater than or equal to the argument key.
 *
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
SpatialKey::operator >= (
   const SpatialKey& spatial_key) const
{
   return ((*this) > spatial_key) || ((*this) == spatial_key);
}

}
}
