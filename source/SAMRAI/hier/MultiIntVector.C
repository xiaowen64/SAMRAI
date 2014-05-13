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
/*
MultiIntVector::MultiIntVector(
   const IntVector& ratio,
   const BlockId& block_id):
   d_vector(1, ratio)
{
   if (s_max_blocks < block_id.getBlockValue() + 1) {
      s_max_blocks = block_id.getBlockValue() + 1;
   }
   for (int b = 0; b < s_max_blocks; ++b) {
      d_vector[b] = hier::IntVector::getOne(ratio.getDim());
   }
   d_vector[block_id.getBlockValue()] = ratio;
}
*/
MultiIntVector::MultiIntVector(
   const IntVector& ratio):
   d_vector(1, ratio)
{
   TBOX_ASSERT(s_max_blocks >= 1);
   if (d_vector.size() != s_max_blocks) {
      d_vector.resize(s_max_blocks, ratio);
   }
   for (int b = 0; b < s_max_blocks; ++b) {
      d_vector[b] = ratio;
   }
}

MultiIntVector::MultiIntVector(
   const std::vector<IntVector>& ratio):
   d_vector(ratio)
{
   if (d_vector.size() > s_max_blocks) {
      s_max_blocks = d_vector.size();
   }
}

MultiIntVector::MultiIntVector(
   const tbox::Dimension& dim,
   int value)
{
   TBOX_ASSERT(s_max_blocks >= 1);
   IntVector tmp(dim, value);
   if (d_vector.size() != s_max_blocks) {
      d_vector.resize(s_max_blocks, tmp);
   }
   for (int b = 0; b < s_max_blocks; ++b) {
      d_vector[b] = tmp;
   }
}

MultiIntVector::MultiIntVector(
   const tbox::Dimension& dim,
   int value,
   int nblocks)
{
   if (nblocks > s_max_blocks) {
      s_max_blocks = nblocks; 
   }
   IntVector tmp(dim, value);
   d_vector.resize(nblocks, tmp);
   for (int b = 0; b < s_max_blocks; ++b) {
      d_vector[b] = tmp;
   }
}



MultiIntVector::MultiIntVector(
   const MultiIntVector& rhs):
   d_vector(rhs.d_vector)
{
}

MultiIntVector::~MultiIntVector()
{
}

std::ostream& operator << (
   std::ostream& s,
   const MultiIntVector& rhs)
{
   s << '(';

   for (int b = 0; b < rhs.d_vector.size(); ++b) {
      s << rhs.d_vector[b];
      if (b < rhs.d_vector.size() - 1)
         s << ",";
   }
   s << ')';

   return s;
}



}
}
