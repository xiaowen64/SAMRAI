/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A box structure representing a portion of the AMR index space
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& DatabaseBox::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
DatabaseBox& DatabaseBox::operator = (
   const DatabaseBox& box)
{
   d_dim = box.d_dim;
   d_data.d_dimension = box.d_data.d_dimension;
   d_data.d_lo[0] = box.d_data.d_lo[0];
   d_data.d_lo[1] = box.d_data.d_lo[1];
   d_data.d_lo[2] = box.d_data.d_lo[2];
   d_data.d_hi[0] = box.d_data.d_hi[0];
   d_data.d_hi[1] = box.d_data.d_hi[1];
   d_data.d_hi[2] = box.d_data.d_hi[2];
   return *this;
}

SAMRAI_INLINE_KEYWORD
void DatabaseBox::setDim(
   const tbox::Dimension& dim)
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   TBOX_ASSERT(dim.getValue() <= DatabaseBox_MAX_DIM);
   d_dim = tbox::Dimension(dim);
   d_data.d_dimension = d_dim.getValue();
}

SAMRAI_INLINE_KEYWORD
int& DatabaseBox::lower(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
   return d_data.d_lo[i];
}

SAMRAI_INLINE_KEYWORD
int& DatabaseBox::upper(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));

   return d_data.d_hi[i];
}

SAMRAI_INLINE_KEYWORD
int DatabaseBox::lower(
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));

   return d_data.d_lo[i];
}

SAMRAI_INLINE_KEYWORD
int DatabaseBox::upper(
   const int i) const
{

   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));

   return d_data.d_hi[i];
}

}
}
