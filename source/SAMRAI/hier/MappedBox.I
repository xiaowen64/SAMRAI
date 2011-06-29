/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Box in the distribued box graph. 
 *
 ************************************************************************/
#include "SAMRAI/hier/MappedBox.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& MappedBox::getDim() const
{
   return d_box.getDim();
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const tbox::Dimension& dim):
   d_id(),
   d_box(dim)
{
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const hier::Box& box,
   const LocalId &local_id,
   const int owner,
   const BlockId &block_id,
   const PeriodicId &periodic_id):
   d_id(local_id, owner, block_id, periodic_id),
   d_box(box)
{
   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() <
      PeriodicShiftCatalog::getCatalog(box.getDim())->getNumberOfShifts());
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const tbox::Dimension& dim,
   const LocalId &local_id,
   const int owner,
   const BlockId &block_id,
   const PeriodicId &periodic_id):
   d_id(local_id, owner, block_id, periodic_id),
   d_box(dim)
{
   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < PeriodicShiftCatalog::getCatalog(
         dim)->getNumberOfShifts());
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const tbox::Dimension& dim,
   const GlobalId& global_id,
   const BlockId &block_id,
   const PeriodicId &periodic_id):
   d_id(global_id, block_id, periodic_id),
   d_box(dim)
{
   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < PeriodicShiftCatalog::getCatalog(
         dim)->getNumberOfShifts());
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const tbox::Dimension& dim,
   const MappedBoxId& mapped_box_id):
   d_id(mapped_box_id),
   d_box(dim)
{
   TBOX_ASSERT(mapped_box_id.getPeriodicId().isValid());
   TBOX_ASSERT(
      mapped_box_id.getPeriodicId().getPeriodicValue() <
      PeriodicShiftCatalog::getCatalog(dim)->getNumberOfShifts());
}

SAMRAI_INLINE_KEYWORD
MappedBox::MappedBox(
   const MappedBox& r):
   d_id(r.d_id),
   d_box(r.d_box)
{
}

SAMRAI_INLINE_KEYWORD
void MappedBox::initialize(
   const hier::Box& box,
   const LocalId &local_id,
   const int owner,
   const BlockId &block_id,
   const PeriodicId &periodic_id)
{
   d_id.initialize(local_id, owner, block_id, periodic_id);
   d_box = box;
}

SAMRAI_INLINE_KEYWORD
hier::MappedBoxId& MappedBox::getId() {
   return d_id;
}

SAMRAI_INLINE_KEYWORD
const hier::MappedBoxId& MappedBox::getId() const {
   return d_id;
}

SAMRAI_INLINE_KEYWORD
hier::Box& MappedBox::getBox() {
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const hier::Box& MappedBox::getBox() const {
   return d_box;
}


SAMRAI_INLINE_KEYWORD
int MappedBox::getOwnerRank() const
{
   return d_id.getOwnerRank();
}


SAMRAI_INLINE_KEYWORD
const LocalId &MappedBox::getLocalId() const
{
   return d_id.getLocalId();
}


SAMRAI_INLINE_KEYWORD
const GlobalId& MappedBox::getGlobalId() const
{
   return d_id.getGlobalId();
}


SAMRAI_INLINE_KEYWORD
const BlockId& MappedBox::getBlockId() const
{
   return d_id.getBlockId();
}


SAMRAI_INLINE_KEYWORD
const PeriodicId &MappedBox::getPeriodicId() const
{
   return d_id.getPeriodicId();
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::isPeriodicImage() const
{
   return d_id.isPeriodicImage();
}

SAMRAI_INLINE_KEYWORD
int MappedBox::commBufferSize(
   const tbox::Dimension& dim)
{
   return MappedBoxId::commBufferSize() + 2 * dim.getValue();
}

SAMRAI_INLINE_KEYWORD
void MappedBox::putToIntBuffer(
   int* buffer) const
{
   d_id.putToIntBuffer(buffer);
   buffer += MappedBoxId::commBufferSize();

   const tbox::Dimension& dim(d_box.getDim());
   for (int d = 0; d < dim.getValue(); ++d) {
      buffer[d * 2] = d_box.lower() (d);
      buffer[d * 2 + 1] = d_box.upper() (d);
   }
}

SAMRAI_INLINE_KEYWORD
void MappedBox::getFromIntBuffer(
   const int* buffer)
{
   d_id.getFromIntBuffer(buffer);
   buffer += MappedBoxId::commBufferSize();

   const tbox::Dimension& dim(d_box.getDim());
   for (int d = 0; d < dim.getValue(); ++d) {
      d_box.lower() (d) = buffer[d * 2];
      d_box.upper() (d) = buffer[d * 2 + 1];
   }
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator == ( const MappedBox& r) const
{
   return d_id == r.d_id;
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator != ( const MappedBox& r) const
{
   return d_id != r.d_id;
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator < ( const MappedBox& r) const
{
   return d_id < r.d_id;
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator > ( const MappedBox& r) const
{
   return d_id > r.d_id;
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator <= ( const MappedBox& r) const
{
   return d_id <= r.d_id;
}


SAMRAI_INLINE_KEYWORD
bool MappedBox::operator >= ( const MappedBox& r) const
{
   return d_id >= r.d_id;
}

}
}
