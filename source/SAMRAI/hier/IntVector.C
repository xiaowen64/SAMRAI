/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/

#ifndef included_hier_IntVector_C
#define included_hier_IntVector_C

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/IntVector.I"
#endif

namespace SAMRAI {
namespace hier {

IntVector * IntVector::s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
IntVector * IntVector::s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

tbox::StartupShutdownManager::Handler
IntVector::s_initialize_finalize_handler(
   IntVector::initializeCallback,
   0,
   0,
   IntVector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

std::istream& operator >> (
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

void IntVector::putUnregisteredToDatabase(
   tbox::Database& database,
   const std::string& name) const
{
   database.putIntegerArray(name, d_vector, d_dim.getValue());
}

void IntVector::getFromDatabase(
   tbox::Database& database,
   const std::string& name)
{
   int d = database.getArraySize(name);
   d_dim = tbox::Dimension(static_cast<unsigned short>(d));
   database.getIntegerArray(name, d_vector, d_dim.getValue());
}

void IntVector::initializeCallback()
{
   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_zeros[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
   }

   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_ones[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);
   }
}

void IntVector::finalizeCallback()
{
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      delete s_zeros[d];
      delete s_ones[d];
   }
}

}
}

#endif
