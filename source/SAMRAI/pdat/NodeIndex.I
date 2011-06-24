/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
void NodeIndex::setOffsets()
{
   const tbox::Dimension& dim(getDim());
   int dim_index = dim.getValue() - 1;
   if (!s_offsets_are_set[dim_index]) {
      s_offsets[dim_index] = std::vector<hier::IntVector>(
            2 << tbox::Dimension::MAXIMUM_DIMENSION_VALUE,
            hier::IntVector(dim));
      for (int i = 0; i < (1 << dim.getValue()); i++) {
         hier::IntVector offset(dim, 0);

         offset(0) = i % 2;
         for (int j = 1; j < dim.getValue(); j++) {
            offset(j) = (i / (1 << j)) % 2;
         }
         s_offsets[dim_index][i] = offset;
      }
      s_offsets_are_set[dim_index] = true;
   }
}

SAMRAI_INLINE_KEYWORD
NodeIndex::NodeIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
   setOffsets();
}

SAMRAI_INLINE_KEYWORD
NodeIndex::NodeIndex(
   const hier::Index& rhs,
   const Corner corner):
   hier::Index(rhs.getDim())
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   setOffsets();
   hier::IntVector::operator = (
      rhs + s_offsets[getDim().getValue() - 1][(int)corner]);
}

SAMRAI_INLINE_KEYWORD
NodeIndex::NodeIndex(
   const hier::Index& rhs,
   const hier::IntVector& corner):
   hier::Index(rhs.getDim())
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < getDim().getValue(); i++) {
      TBOX_ASSERT(corner(i) == 0 || corner(i) == 1);
   }
#endif
   setOffsets();
   hier::IntVector::operator = (
      rhs + corner);
}

SAMRAI_INLINE_KEYWORD
NodeIndex::NodeIndex(
   const NodeIndex& rhs):
   hier::Index(rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   setOffsets();
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator = (
   const NodeIndex& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator = (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex::~NodeIndex()
{
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   NodeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator + (
   const int rhs) const
{
   NodeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   NodeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator - (
   const int rhs) const
{
   NodeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   NodeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
NodeIndex& NodeIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
NodeIndex NodeIndex::operator * (
   const int rhs) const
{
   NodeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool NodeIndex::operator == (
   const NodeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator == (rhs);
}

SAMRAI_INLINE_KEYWORD
bool NodeIndex::operator != (
   const NodeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator != (rhs);
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
