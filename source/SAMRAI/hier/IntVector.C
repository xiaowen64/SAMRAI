/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

IntVector * IntVector::s_zeros[SAMRAI::MAX_DIM_VAL];
IntVector * IntVector::s_ones[SAMRAI::MAX_DIM_VAL];

tbox::StartupShutdownManager::Handler
IntVector::s_initialize_finalize_handler(
   IntVector::initializeCallback,
   0,
   0,
   IntVector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

IntVector::IntVector(
   const tbox::Dimension& dim):
   d_dim(dim)
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif

}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int value):
   d_dim(dim)
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = value;
   }

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = d_dim.getValue(); i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const std::vector<int>& a):
   d_dim(static_cast<unsigned short>(a.size()))
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = a[i];
   }

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = d_dim.getValue(); i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const IntVector& rhs):
   d_dim(rhs.getDim())
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = rhs.d_vector[i];
   }
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int array[]):
   d_dim(dim)
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = array[i];
   }
}

IntVector::~IntVector()
{
}

std::istream&
operator >> (
   std::istream& s,
   IntVector& rhs)
{
   while (s.get() != '(') ;

   for (int i = 0; i < rhs.getDim().getValue(); i++) {
      s >> rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         while (s.get() != ',') ;
   }

   while (s.get() != ')') ;

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const IntVector& rhs)
{
   s << '(';

   for (int i = 0; i < rhs.getDim().getValue(); i++) {
      s << rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         s << ",";
   }
   s << ')';

   return s;
}

void
IntVector::putToRestart(
   tbox::Database& restart_db,
   const std::string& name) const
{
   restart_db.putIntegerArray(name, d_vector, d_dim.getValue());
}

void
IntVector::getFromRestart(
   tbox::Database& restart_db,
   const std::string& name)
{
   TBOX_ASSERT(d_dim.getValue() ==
      static_cast<unsigned short>(restart_db.getArraySize(name)));
   restart_db.getIntegerArray(name, d_vector, d_dim.getValue());
}

/*
 *************************************************************************
 * Sort the sizes of the given IntVector from smallest to largest value.
 *************************************************************************
 */
void
IntVector::sortIntVector(
   const IntVector& values)
{
   const IntVector num_cells = values;

   for (int d = 0; d < d_dim.getValue(); d++) {
      d_vector[d] = d;
   }
   for (int d0 = 0; d0 < d_dim.getValue() - 1; d0++) {
      for (int d1 = d0 + 1; d1 < d_dim.getValue(); d1++) {
         if (values(d_vector[d0]) > values(d_vector[d1])) {
            int tmp_d = d_vector[d0];
            d_vector[d0] = d_vector[d1];
            d_vector[d1] = tmp_d;
         }
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int d = 0; d < d_dim.getValue() - 1; d++) {
      TBOX_ASSERT(values(d_vector[d]) <= values(d_vector[d + 1]));
   }
#endif
}

void
IntVector::initializeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_zeros[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
      s_ones[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);
   }
}

void
IntVector::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      delete s_zeros[d];
      delete s_ones[d];
   }
}

}
}
