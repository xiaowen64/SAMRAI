/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   A container of IntVectors
 *
 ************************************************************************/
#include "SAMRAI/hier/MultiIntVector.h"
#if 0
namespace SAMRAI {
namespace hier {

int MultiIntVector::s_num_blocks = 1;


MultiIntVector::MultiIntVector(
   const std::vector<IntVector>& vector):
   d_vector(vector),
   d_dim(vector[0].getDim())
{
   TBOX_ASSERT(!vector.empty());
   TBOX_ASSERT(s_num_blocks == 1 || s_num_blocks == d_vector.size());
   if (s_num_blocks == 1) {
      s_num_blocks = d_vector.size();
   }
}

MultiIntVector::MultiIntVector(
   const IntVector& vector):
   d_vector(1, vector),
   d_dim(vector.getDim())
{
   TBOX_ASSERT(s_num_blocks > 0);
   if (d_vector.size() != s_num_blocks) {
      d_vector.resize(s_num_blocks, vector);
   }
   for (int b = 0; b < s_num_blocks; ++b) {
      d_vector[b] = vector;
   }
}

MultiIntVector::MultiIntVector(
   const tbox::Dimension& dim,
   int value):
   d_vector(1, IntVector(dim, value)),
   d_dim(dim)
{
   TBOX_ASSERT(s_num_blocks > 0);
   if (d_vector.size() != s_num_blocks) {
      d_vector.resize(s_num_blocks, d_vector[0]);
   }
   for (int b = 1; b < s_num_blocks; ++b) {
      d_vector[b] = d_vector[0];
   }
}

MultiIntVector::MultiIntVector(
   const MultiIntVector& copy_vector):
   d_vector(copy_vector.d_vector),
   d_dim(copy_vector.d_dim)
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

std::istream&
operator >> (
   std::istream& s,
   MultiIntVector& rhs)
{
   while (s.get() != '(') ;

   for (int b = 0; b < rhs.d_vector.size(); ++b) {
      s >> rhs.d_vector[b];
      if (b < rhs.d_vector.size() - 1)
         while (s.get() != ',') ;
   }

   while (s.get() != ')') ;

   return s;
}

}
}
#endif
