/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility for computing all possible periodic shift directions.
 *
 ************************************************************************/
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const PeriodicId& PeriodicShiftCatalog::getZeroShiftNumber() const
{
   return PeriodicId::zero();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PeriodicShiftCatalog::isPeriodic() const
{
   return d_shifts.size() > 1;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
int PeriodicShiftCatalog::getNumberOfShifts() const
{
   return (int)d_shifts.size();
}

SAMRAI_INLINE_KEYWORD
const IntVector& PeriodicShiftCatalog::shiftNumberToShiftDistance(
   const PeriodicId& shift_number) const
{
   TBOX_ASSERT(shift_number.isValid());
   TBOX_ASSERT(shift_number.getPeriodicValue() < (int)d_shifts.size());

   return d_shifts[shift_number.getPeriodicValue()];
}

SAMRAI_INLINE_KEYWORD
const PeriodicId& PeriodicShiftCatalog::getOppositeShiftNumber(
   const PeriodicId& shift_number) const
{
   TBOX_ASSERT(shift_number.isValid());
   TBOX_ASSERT(shift_number.getPeriodicValue() < (int)d_shifts.size());

   return d_opposite_number[shift_number.getPeriodicValue()];
}

SAMRAI_INLINE_KEYWORD
PeriodicId PeriodicShiftCatalog::shiftDistanceToShiftNumber(
   const IntVector& shift_distance) const
{
   unsigned int s;
   for (s = 0; s < d_shifts.size(); ++s) {
      if (d_shifts[s] == shift_distance) break;
   }
   if (s == d_shifts.size()) {
      s = tbox::MathUtilities<int>::getMax();
   }
   return PeriodicId(s);
}

SAMRAI_INLINE_KEYWORD
int PeriodicShiftCatalog::getInvalidShiftNumber() const
{
   return tbox::MathUtilities<int>::getMax();
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
