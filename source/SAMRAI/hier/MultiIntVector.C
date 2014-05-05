/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/
#include "SAMRAI/hier/MultiIntVector.h"

namespace SAMRAI {
namespace hier {

int MultiIntVector::s_max_blocks = 0;

MultiIntVector::MultiIntVector(
   const IntVector& ratio,
   const BlockId& block_id):
   d_ratio(1, ratio),
   d_is_multiblock(false)
{
   if (s_max_blocks < block_id.getBlockValue() + 1) {
      s_max_blocks = block_id.getBlockValue() + 1;
   }
   for (int b = 0; b < s_max_blocks; ++b) {
      d_ratio[b] = hier::IntVector::getOne(ratio.getDim());
   }
   d_ratio[block_id.getBlockValue()] = ratio;
}

MultiIntVector::MultiIntVector(
   const IntVector& ratio):
   d_ratio(1, ratio),
   d_is_multiblock(false)
{
   TBOX_ASSERT(s_max_blocks >= 1);
   for (int b = 0; b < s_max_blocks; ++b) {
      d_ratio[b] = ratio;
   }
}

MultiIntVector::MultiIntVector(
   const std::vector<IntVector>& ratio):
   d_ratio(ratio),
   d_is_multiblock(false)
{
   if (d_ratio.size() > s_max_blocks) {
      s_max_blocks = d_ratio.size();
   }
   if (d_ratio.size() > 1) {
      d_is_multiblock = true;
   }
}

MultiIntVector::MultiIntVector(
   const MultiIntVector& rhs):
   d_ratio(rhs.d_ratio),
   d_is_multiblock(rhs.d_is_multiblock)
{
}

MultiIntVector::~MultiIntVector()
{
}


}
}
