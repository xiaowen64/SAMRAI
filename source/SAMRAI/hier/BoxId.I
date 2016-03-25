/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/

#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxId::BoxId():
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
BoxId::BoxId(
   const GlobalId& id,
   const BlockId& block_id,
   const PeriodicId& periodic_id):
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
BoxId::BoxId(
   const LocalId& local_id,
   const int owner,
   const BlockId& block_id,
   const PeriodicId& periodic_id):
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
BoxId::BoxId(
   const BoxId& r):
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
BoxId::~BoxId()
{
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxId::initialize(
   const LocalId& local_id,
   const int owner,
   const BlockId& block_id,
   const PeriodicId& periodic_id)
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
const hier::GlobalId& BoxId::getGlobalId() const {
   return d_global_id;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int BoxId::getOwnerRank() const {
   return d_global_id.getOwnerRank();
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const LocalId& BoxId::getLocalId() const {
   return d_global_id.getLocalId();
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const BlockId& BoxId::getBlockId() const {
   return d_block_id;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const PeriodicId& BoxId::getPeriodicId() const {
   return d_periodic_id;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::isPeriodicImage() const {
   return d_periodic_id != PeriodicId::zero();
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int BoxId::commBufferSize()
{
   return 4;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxId::putToIntBuffer(
   int* buffer) const
{
   buffer[0] = getOwnerRank();
   buffer[1] = getLocalId().getValue();
   buffer[2] = getBlockId().getBlockValue();
   buffer[3] = getPeriodicId().getPeriodicValue();
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxId::getFromIntBuffer(
   const int* buffer)
{
   initialize(LocalId(buffer[1]),
      buffer[0],
      BlockId(buffer[2]),
      PeriodicId(buffer[3]));
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::operator == (
   const BoxId& r) const {
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
bool BoxId::operator != (
   const BoxId& r) const {
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
bool BoxId::operator < (
   const BoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return d_global_id.getOwnerRank() < r.d_global_id.getOwnerRank() ||
          (d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
           (d_block_id < r.d_block_id ||
            (d_block_id == r.d_block_id &&
             (d_global_id.getLocalId() < r.d_global_id.getLocalId() ||
              (d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
               (d_periodic_id < r.d_periodic_id))))));
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::operator > (
   const BoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return d_global_id.getOwnerRank() > r.d_global_id.getOwnerRank() ||
          (d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
           (d_block_id > r.d_block_id ||
            (d_block_id == r.d_block_id &&
             (d_global_id.getLocalId() > r.d_global_id.getLocalId() ||
              (d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
               (d_periodic_id > r.d_periodic_id))))));
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::operator <= (
   const BoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return *this < r || *this == r;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::operator >= (
   const BoxId& r) const {
   TBOX_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
   return *this > r || *this == r;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxId::isValid() const
{
   return (d_block_id.isValid() &&
           d_periodic_id.isValid() &&
           d_global_id.getLocalId() != LocalId::getInvalidId() && 
           d_global_id.getOwnerRank() != tbox::SAMRAI_MPI::getInvalidRank());
            
}

}
}
