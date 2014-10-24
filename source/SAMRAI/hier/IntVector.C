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


/*
 * *************************************************************************
 * Constructors
 * *************************************************************************
 */

IntVector::IntVector(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_block_size(1),
   d_vector(dim.getValue())
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   int block_size,
   const tbox::Dimension& dim):
   d_dim(dim),
   d_block_size(block_size),
   d_vector(dim.getValue()*block_size)
{
   TBOX_ASSERT(block_size >=1);
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < block_size*dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   int value,
   int block_size):
   d_dim(dim),
   d_block_size(block_size),
   d_vector(dim.getValue()*block_size, value)
{
   TBOX_ASSERT(block_size >=1);
}

IntVector::IntVector(
   const std::vector<int>& vec,
   int block_size):
   d_dim(static_cast<unsigned short>(vec.size())),
   d_block_size(block_size),
   d_vector(vec.size()*block_size)
{
   TBOX_ASSERT(vec.size() >= 1);
   for (int b = 0; b < block_size; ++b) {
      int offset = b*d_dim.getValue();
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = vec[i];
      }
   }
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int array[],
   int block_size):
   d_dim(dim),
   d_block_size(block_size),
   d_vector(dim.getValue()*block_size)
{
   for (int b = 0; b < block_size; ++b) {
      int offset = b*d_dim.getValue();
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = array[i];
      }
   }
}

IntVector::IntVector(
   const IntVector& rhs):
   d_dim(rhs.getDim()),
   d_block_size(rhs.d_block_size),
   d_vector(rhs.d_vector)
{
   TBOX_ASSERT(d_block_size >= 1);
}

IntVector::IntVector(
   const IntVector& rhs,
   int block_size):
   d_dim(rhs.getDim()),
   d_block_size(block_size),
   d_vector(rhs.getDim().getValue() * block_size)
{
   TBOX_ASSERT(d_block_size >= 1);
   TBOX_ASSERT(rhs.d_block_size == d_block_size || rhs.d_block_size == 1); 
   if (rhs.d_block_size == 1 && d_block_size != 1) {
      for (int b = 0; b < d_block_size; ++b) {
         int offset = b*d_dim.getValue();
         for (int i = 0; i < d_dim.getValue(); ++i) {
            d_vector[offset + i] = rhs.d_vector[i];
         }
      }
   } else {
      d_vector = rhs.d_vector;
   }
}

IntVector::IntVector(
   const Index& rhs,
   int block_size):
   d_dim(rhs.getDim()),
   d_block_size(block_size),
   d_vector(rhs.getDim().getValue() * block_size)
{
   TBOX_ASSERT(d_block_size >= 1);
   for (int b = 0; b < block_size; ++b) {
      int offset = b*d_dim.getValue();
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = rhs[i];
      } 
   }
}

/*
 * *************************************************************************
 * Destructor 
 * *************************************************************************
 */
IntVector::~IntVector()
{
}

/*
 * *************************************************************************
 * Assignment
 * *************************************************************************
 */
IntVector&
IntVector::operator = (
   const Index& rhs)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
   if (d_block_size != 1) {
      d_block_size = 1;
      d_vector.resize(d_dim.getValue());
   }

   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = rhs[i];
   }
   return *this;
}

/*
 * *************************************************************************
 * Streaming I/O
 * *************************************************************************
 */
std::istream&
operator >> (
   std::istream& s,
   IntVector& rhs)
{
   for (int b = 0; b < rhs.getBlockSize(); ++b) {
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

   for (int b = 0; b < rhs.getBlockSize(); ++b) {
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

/*
 * *************************************************************************
 * Write/read for restart
 * *************************************************************************
 */
void
IntVector::putToRestart(
   tbox::Database& restart_db,
   const std::string& name) const
{
   boost::shared_ptr<tbox::Database> intvec_db =
      restart_db.putDatabase(name);
   intvec_db->putInteger("d_block_size", d_block_size);
   intvec_db->putIntegerVector("d_vector",
                               d_vector);

}

void
IntVector::getFromRestart(
   tbox::Database& restart_db,
   const std::string& name)
{
   boost::shared_ptr<tbox::Database> intvec_db =
      restart_db.getDatabase(name);

   d_block_size = intvec_db->getInteger("d_block_size");
   d_vector = intvec_db->getIntegerVector("d_vector"); 

   TBOX_ASSERT(d_block_size * d_dim.getValue() == d_vector.size());

}

/*
 *************************************************************************
 * Sort the values of the given IntVector from smallest to largest value.
 *************************************************************************
 */
void
IntVector::sortIntVector(
   const IntVector& values)
{
   for (int b = 0; b < d_block_size; ++b ) {
      int offset = b*d_dim.getValue();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         d_vector[offset + d] = d;
      }
      for (int d0 = 0; d0 < d_dim.getValue() - 1; ++d0) {
         for (int d1 = d0 + 1; d1 < d_dim.getValue(); ++d1) {
            if (values(d_vector[offset + d0]) > values(d_vector[offset + d1])) {
               int tmp_d = d_vector[offset + d0];
               d_vector[offset + d0] = d_vector[offset + d1];
               d_vector[offset + d1] = tmp_d;
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (int d = 0; d < d_dim.getValue() - 1; ++d) {
         TBOX_ASSERT(values(d_vector[offset + d]) <= values(d_vector[offset + d + 1]));
      }
#endif
   }
}

/*
 * *************************************************************************
 * Callback routines
 * *************************************************************************
 */
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
