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

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"


namespace SAMRAI {
namespace hier {

IntVector * IntVector::s_zeros[SAMRAI::MAX_DIM_VAL];
IntVector * IntVector::s_ones[SAMRAI::MAX_DIM_VAL];
IntVector * IntVector::s_mb_zeros[SAMRAI::MAX_DIM_VAL];
IntVector * IntVector::s_mb_ones[SAMRAI::MAX_DIM_VAL];

tbox::StartupShutdownManager::Handler
IntVector::s_initialize_finalize_handler(
   IntVector::initializeCallback,
   0,
   0,
   IntVector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

IntVector::IntVector(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_size(1),
   d_vector(new int[dim.getValue()])
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   int size,
   const tbox::Dimension& dim):
   d_dim(dim),
   d_size(size),
   d_vector(new int[dim.getValue()*size])
{
   TBOX_ASSERT(size >=1);
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < size*dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   int value,
   int size):
   d_dim(dim),
   d_size(size),
   d_vector(new int[dim.getValue()*size])
{
   TBOX_ASSERT(size >=1);
   for (int i = 0; i < size*dim.getValue(); ++i) {
      d_vector[i] = value;
   }
}

IntVector::IntVector(
   const std::vector<int>& a,
   int size):
   d_dim(static_cast<unsigned short>(a.size())),
   d_size(size),
   d_vector(new int[a.size()*size])
{
   for (int b = 0; b < size; ++b) {
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[b*d_dim.getValue() + i] = a[i];
      }
   }
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int array[],
   int size):
   d_dim(dim),
   d_size(size),
   d_vector(new int[dim.getValue()*size])
{
   for (int b = 0; b < size; ++b) {
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[b*d_dim.getValue() + i] = array[i];
      }
   }
}


IntVector::IntVector(
   const IntVector& rhs):
   d_dim(rhs.getDim()),
   d_size(rhs.d_size),
   d_vector(new int[rhs.getDim().getValue() * rhs.d_size])
{
   TBOX_ASSERT(d_size >= 1);
   for (int b = 0; b < d_size; ++b) {
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[b*d_dim.getValue()+i] = rhs.d_vector[b*d_dim.getValue()+i];
      }
   }
}

IntVector::IntVector(
   const IntVector& rhs,
   int size):
   d_dim(rhs.getDim()),
   d_size(size),
   d_vector(new int[rhs.getDim().getValue() * size])
{
   TBOX_ASSERT(d_size >= 1);
   TBOX_ASSERT(rhs.d_size == d_size || rhs.d_size == 1); 
   if (rhs.d_size == 1) {
      for (int b = 0; b < d_size; ++b) {
         for (int i = 0; i < d_dim.getValue(); ++i) {
            d_vector[b*d_dim.getValue() + i] = rhs.d_vector[i];
         }
      }
   } else {
      for (int b = 0; b < d_size; ++b) {
         for (int i = 0; i < d_dim.getValue(); ++i) {
            d_vector[b*d_dim.getValue()+i] = rhs.d_vector[b*d_dim.getValue()+i];
         }
      }
   }
}


IntVector::IntVector(
   const Index& rhs,
   int size):
   d_dim(rhs.getDim()),
   d_size(size),
   d_vector(new int[rhs.getDim().getValue() * size])
{
   TBOX_ASSERT(d_size >= 1);
   for (int b = 0; b < size; ++b) {
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[b*d_dim.getValue() + i] = rhs[i];
      }
   }
}

IntVector::~IntVector()
{
   if (d_vector) {
      delete[] d_vector;
   }
}

IntVector&
IntVector::operator = (
   const Index& rhs)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
   TBOX_ASSERT(d_size >= 1);
   if (d_size != 1) {
      delete[] d_vector;
      d_size = 1;
      d_vector = new int[d_dim.getValue()];
   }

   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = rhs[i];
   }
   return *this;
}

std::istream&
operator >> (
   std::istream& s,
   IntVector& rhs)
{
   for (int b = 0; b < rhs.size(); ++b) {
      while (s.get() != '(') ;
      for (int i = 0; i < rhs.getDim().getValue(); ++i) {
         s >> rhs(b,i);
         if (i < rhs.getDim().getValue() - 1)
            while (s.get() != ',') ;
      }
      while (s.get() != ')') ;
   }

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const IntVector& rhs)
{

   for (int b = 0; b < rhs.size(); ++b) {
      s << '(';
      for (int i = 0; i < rhs.getDim().getValue(); ++i) {
         s << rhs(b,i);
         if (i < rhs.getDim().getValue() - 1)
            s << ",";
      }
      s << ')';
   }

   return s;
}

void
IntVector::putToRestart(
   tbox::Database& restart_db,
   const std::string& name) const
{
   boost::shared_ptr<tbox::Database> intvec_db =
      restart_db.putDatabase(name);
   intvec_db->putInteger("d_size", d_size);
   intvec_db->putIntegerArray("d_vector",
                              d_vector,
                              d_size * d_dim.getValue());

}

void
IntVector::getFromRestart(
   tbox::Database& restart_db,
   const std::string& name)
{
   boost::shared_ptr<tbox::Database> intvec_db =
      restart_db.getDatabase(name);

   int size = intvec_db->getInteger("d_size");
   std::vector<int> tmp_vector(size * d_dim.getValue());
   TBOX_ASSERT(size * d_dim.getValue() == intvec_db->getArraySize("d_vector"));

   if (size != d_size) {
      delete[] d_vector;
      d_size = size;
      d_vector = new int[d_dim.getValue()*d_size];
   }

   intvec_db->getIntegerArray("d_vector",
                              d_vector,
                              d_size * d_dim.getValue());

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
   for (int b = 0; b < d_size; ++b ) {
      for (int d = 0; d < d_dim.getValue(); ++d) {
         d_vector[b*d_dim.getValue() + d] = d;
      }
      for (int d0 = 0; d0 < d_dim.getValue() - 1; ++d0) {
         for (int d1 = d0 + 1; d1 < d_dim.getValue(); ++d1) {
            if (values(d_vector[b*d_dim.getValue() + d0]) > values(d_vector[b*d_dim.getValue() + d1])) {
               int tmp_d = d_vector[b*d_dim.getValue() + d0];
               d_vector[b*d_dim.getValue() + d0] = d_vector[b*d_dim.getValue() + d1];
               d_vector[b*d_dim.getValue() + d1] = tmp_d;
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (int d = 0; d < d_dim.getValue() - 1; ++d) {
         TBOX_ASSERT(values(d_vector[b*d_dim.getValue() + d]) <= values(d_vector[b*d_dim.getValue() + d + 1]));
      }
#endif
   }
}

void
IntVector::initializeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_zeros[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
      s_ones[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);
      s_mb_zeros[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
      s_mb_ones[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);
   }
}

void
IntVector::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      delete s_zeros[d];
      delete s_ones[d];
      delete s_mb_zeros[d];
      delete s_mb_ones[d];
   }
}

}
}
