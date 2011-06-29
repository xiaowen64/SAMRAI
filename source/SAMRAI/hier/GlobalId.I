/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Globally unique identifier that can be locally determined.
 *
 ************************************************************************/
#include "SAMRAI/hier/GlobalId.h"

#ifndef included_tbox_Utilities
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace hier {

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
GlobalId::GlobalId():
   d_owner_rank(tbox::SAMRAI_MPI::getInvalidRank()),
   d_local_id(LocalId::getInvalidId())
{
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
GlobalId::GlobalId(
   const LocalId &local_id,
   const int owner):
   d_owner_rank(owner),
   d_local_id(local_id)
{
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
GlobalId::GlobalId(
   const GlobalId& other):
   d_owner_rank(other.d_owner_rank),
   d_local_id(other.d_local_id)
{
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
GlobalId::~GlobalId()
{
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
int& GlobalId::getOwnerRank() {
   return d_owner_rank;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const int& GlobalId::getOwnerRank() const {
   return d_owner_rank;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId& GlobalId::getLocalId() {
   return d_local_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const LocalId& GlobalId::getLocalId() const {
   return d_local_id;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator == (
   const GlobalId& r) const {
   bool rval = (d_owner_rank == r.d_owner_rank) &&
      (d_local_id == r.d_local_id);
   return rval;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator != (
   const GlobalId& r) const {
   bool rval = (d_owner_rank != r.d_owner_rank) ||
      (d_local_id != r.d_local_id);
   return rval;
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator < (
   const GlobalId& r) const {
   return (d_owner_rank < r.d_owner_rank) ||
          ((d_owner_rank == r.d_owner_rank) && (d_local_id < r.d_local_id));
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator > (
   const GlobalId& r) const {
   return (d_owner_rank > r.d_owner_rank) ||
          ((d_owner_rank == r.d_owner_rank) && (d_local_id > r.d_local_id));
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator <= (
   const GlobalId& r) const {
   return (d_owner_rank < r.d_owner_rank) ||
          ((d_owner_rank == r.d_owner_rank) &&
           (d_local_id <= r.d_local_id));
}

/*
************************************************************************
************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool GlobalId::operator >= (
   const GlobalId& r) const {
   return (d_owner_rank > r.d_owner_rank) ||
          ((d_owner_rank == r.d_owner_rank) &&
           (d_local_id >= r.d_local_id));
}

}
}
