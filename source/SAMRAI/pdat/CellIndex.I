/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
CellIndex::CellIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

SAMRAI_INLINE_KEYWORD
CellIndex::CellIndex(
   const hier::Index& rhs):hier::Index(rhs)
{
}

SAMRAI_INLINE_KEYWORD
CellIndex::CellIndex(
   const CellIndex& rhs):hier::Index(rhs)
{
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator = (
   const CellIndex& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   hier::Index::operator = (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex::~CellIndex()
{
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   CellIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator + (
   const int rhs) const
{
   CellIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   CellIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator - (
   const int rhs) const
{
   CellIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   CellIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
CellIndex& CellIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIndex CellIndex::operator * (
   const int rhs) const
{
   CellIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool CellIndex::operator == (
   const CellIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   return ((hier::Index *)this)->operator == (rhs);
}

SAMRAI_INLINE_KEYWORD
bool CellIndex::operator != (
   const CellIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   return ((hier::Index *)this)->operator != (rhs);
}

}
}
