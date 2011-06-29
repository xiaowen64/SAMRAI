/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Box in the distribued box graph. 
 *
 ************************************************************************/
#ifndef included_hier_MappedBox_C
#define included_hier_MappedBox_C

#include "SAMRAI/hier/MappedBox.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"

#include <iostream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBox.I"
#endif

namespace SAMRAI {
namespace hier {

/*
**********************************************************************************
**********************************************************************************
*/
MappedBox::~MappedBox()
{
}


/*
 *********************************************************************************
 * Construct MappedBox from the components of a reference MappedBox
 * and possibly changing the periodic shift.
 *
 * This method is not inlined because constructing a periodic-shifted
 * MappedBox from another (possibly shifted) MappedBox is more involved
 * and less frequently used.
 *********************************************************************************
 */
MappedBox::MappedBox(
   const MappedBox& r,
   const PeriodicId &periodic_id,
   const IntVector& refinement_ratio):
   d_id(r.getLocalId(), r.getOwnerRank(), r.getBlockId(), periodic_id),
   d_box(r.d_box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(d_box, r, refinement_ratio);

   const tbox::Dimension& dim(d_box.getDim());

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < shift_catalog->getNumberOfShifts());

   if (refinement_ratio > hier::IntVector::getZero(dim)) {

      if (r.getPeriodicId() != shift_catalog->getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         d_box.shift(-shift_catalog->shiftNumberToShiftDistance(r.
               getPeriodicId())
            * refinement_ratio);
      }

      if (periodic_id != shift_catalog->getZeroShiftNumber()) {
         // Apply the shift for this MappedBox.
         d_box.shift(shift_catalog->shiftNumberToShiftDistance(periodic_id)
            * refinement_ratio);
      }

   } else if (refinement_ratio < hier::IntVector::getZero(dim)) {

      if (r.getPeriodicId() != shift_catalog->getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         d_box.shift(shift_catalog->shiftNumberToShiftDistance(r.getPeriodicId())
            / refinement_ratio);
      }

      if (periodic_id != shift_catalog->getZeroShiftNumber()) {
         // Apply the shift for this MappedBox.
         d_box.shift(-shift_catalog->shiftNumberToShiftDistance(periodic_id)
            / refinement_ratio);
      }

   } else {

      TBOX_ERROR(
         "MappedBox::MappedBox: Invalid refinement ratio "
         << refinement_ratio
         <<
         "\nRefinement ratio must be completely positive or negative.");

   }
}


/*
**********************************************************************************
* Default constructor needed by STL classes but should not be used
* by anything.
**********************************************************************************
*/
MappedBox::MappedBox():
   d_box(tbox::Dimension::getInvalidDimension())
{
}


/*
 *********************************************************************************
 * Construct MappedBox from a reference MappedBox and possibly
 * changing the periodic shift.
 *
 * This method is not inlined because initializing a periodic-shifted
 * MappedBox from another (possibly shifted) MappedBox is more involved
 * and less frequently used.
 *
 * We inititalize d_id last so that we can support inititalizing an
 * object from a reference to itself.
 *********************************************************************************
 */
void MappedBox::initialize(
   const MappedBox& r,
   const PeriodicId &periodic_id,
   const IntVector& refinement_ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(d_box, r, refinement_ratio);

   const tbox::Dimension& dim(d_box.getDim());

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < shift_catalog->getNumberOfShifts());

   d_box = r.d_box;

   if (refinement_ratio > hier::IntVector::getZero(dim)) {

      if (r.getPeriodicId() != shift_catalog->getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         d_box.shift(-shift_catalog->shiftNumberToShiftDistance(r.
               getPeriodicId())
            * refinement_ratio);
      }

      if (periodic_id != shift_catalog->getZeroShiftNumber()) {
         // Apply the shift for this MappedBox.
         d_box.shift(shift_catalog->shiftNumberToShiftDistance(periodic_id)
            * refinement_ratio);
      }

   } else if (refinement_ratio < hier::IntVector::getZero(dim)) {

      if (r.getPeriodicId() != shift_catalog->getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         d_box.shift(shift_catalog->shiftNumberToShiftDistance(r.getPeriodicId())
            / refinement_ratio);
      }

      if (periodic_id != shift_catalog->getZeroShiftNumber()) {
         // Apply the shift for this MappedBox.
         d_box.shift(-shift_catalog->shiftNumberToShiftDistance(periodic_id)
            / refinement_ratio);
      }

   } else {

      TBOX_ERROR(
         "MappedBox::initialize: Invalid refinement ratio "
         << refinement_ratio
         <<
         "\nRefinement ratio must be completely positive or negative.");

   }

   d_id.initialize(r.getLocalId(), r.getOwnerRank(), r.getBlockId(), periodic_id);
}


/*
*******************************************************************************
Stream-insert operator.
*******************************************************************************
*/
std::ostream& operator << (
   std::ostream& co,
   const MappedBox& r)
{
   co << r.d_id << ':' << r.getBox();
   return co;
}

}
}
#endif
