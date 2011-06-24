/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a MappedBox.
 *
 ************************************************************************/

#include "SAMRAI/hier/MappedBoxId.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
MappedBoxId::MappedBoxId():
   d_global_id(),
   d_block_id(),
   d_periodic_id()
{
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
MappedBoxId::MappedBoxId(
   const GlobalId& id,
   const BlockId &block_id,
   const PeriodicId &periodic_id):
   d_global_id(id),
   d_block_id(block_id),
   d_periodic_id(periodic_id)
{
   TBOX_ASSERT(periodic_id.isValid());
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
MappedBoxId::MappedBoxId(
   const LocalId &local_id,
   const int owner,
   const BlockId &block_id,
   const PeriodicId &periodic_id):
   d_global_id(local_id, owner),
   d_block_id(block_id),
   d_periodic_id(periodic_id)
{
   TBOX_ASSERT(periodic_id.isValid());
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
MappedBoxId::MappedBoxId(
   const MappedBoxId& r):
   d_global_id(r.d_global_id),
   d_block_id(r.d_block_id),
   d_periodic_id(r.d_periodic_id)
{
   TBOX_ASSERT(r.d_periodic_id.isValid());
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
MappedBoxId::~MappedBoxId()
{
   return;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
void MappedBoxId::initialize(
   const LocalId &local_id,
   const int owner,
   const BlockId &block_id,
   const PeriodicId &periodic_id)
{
   TBOX_ASSERT(periodic_id.isValid());

   d_global_id.getLocalId() = local_id;
   d_global_id.getOwnerRank() = owner;
   d_block_id = block_id;
   d_periodic_id = periodic_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const hier::GlobalId& MappedBoxId::getGlobalId() const {
   return d_global_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
int MappedBoxId::getOwnerRank() const {
   return d_global_id.getOwnerRank();
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const LocalId &MappedBoxId::getLocalId() const {
   return d_global_id.getLocalId();
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const BlockId &MappedBoxId::getBlockId() const {
   return d_block_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const PeriodicId &MappedBoxId::getPeriodicId() const {
   return d_periodic_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::isPeriodicImage() const {
   return d_periodic_id != PeriodicId::zero();
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
int MappedBoxId::commBufferSize()
{
   return 4;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
void MappedBoxId::putToIntBuffer(
   int* buffer) const
{
   buffer[0] = getOwnerRank();
   buffer[1] = getLocalId().getValue();
   buffer[2] = getBlockId().getBlockValue();
   buffer[3] = getPeriodicId().getPeriodicValue();
   return;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
void MappedBoxId::getFromIntBuffer(
   const int* buffer)
{
   initialize(LocalId(buffer[1]),
              buffer[0],
              BlockId(buffer[2]),
              PeriodicId(buffer[3]));
   return;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator == (
   const MappedBoxId& r) const {
   bool rval =
      d_block_id == r.d_block_id &&
      d_global_id == r.d_global_id &&
      d_periodic_id == r.d_periodic_id;
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return rval;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator != (
   const MappedBoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   bool rval =
      d_block_id != r.d_block_id ||
      d_global_id != r.d_global_id ||
      d_periodic_id != r.d_periodic_id;
   return rval;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator < (
   const MappedBoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return
      d_global_id.getOwnerRank() < r.d_global_id.getOwnerRank() ||
      ( d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
        ( d_block_id < r.d_block_id ||
          ( d_block_id == r.d_block_id &&
            ( d_global_id.getLocalId() < r.d_global_id.getLocalId() ||
              ( d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
                ( d_periodic_id < r.d_periodic_id ) ) ) ) ) );
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator > (
   const MappedBoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return
      d_global_id.getOwnerRank() > r.d_global_id.getOwnerRank() ||
      ( d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
        ( d_block_id > r.d_block_id ||
          ( d_block_id == r.d_block_id &&
            ( d_global_id.getLocalId() > r.d_global_id.getLocalId() ||
              ( d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
                ( d_periodic_id > r.d_periodic_id ) ) ) ) ) );
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator <= (
   const MappedBoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return *this < r || *this == r;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool MappedBoxId::operator >= (
   const MappedBoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return *this > r || *this == r;
}

}
}
