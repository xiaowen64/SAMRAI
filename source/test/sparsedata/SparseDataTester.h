/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Test class for SparseData.
 *
 ************************************************************************/
#ifndef included_SparseDataTester_h
#define included_SparseDataTester_h

#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/pdat/SparseData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/hier/Index.h"

#include <string>
#include <vector>

namespace sam_test {

using namespace SAMRAI;

class SparseDataTester
{
public:
   SparseDataTester(const SAMRAI::tbox::Dimension& dim);
   ~SparseDataTester();

   bool testConstruction();
   bool testAdd();
   bool testRemove();
   bool testCopy();
   bool testCopy2();
   void testTiming();
   bool testPackStream();
   bool testDatabaseInterface();

private:
   static const int DSIZE = 7;
   static const int ISIZE = 3;

   typedef pdat::SparseData<pdat::CellGeometry> SparseDataType;

   void _fillObject(tbox::Pointer<SparseDataType> sparse_data);
   void _getDblKeys(std::vector<std::string>& keys);
   void _getIntKeys(std::vector<std::string>& keys);
   void _getDblValues(double* values);
   void _getIntValues(int* values);
   bool _testCopy(tbox::Pointer<SparseDataType> src,
      tbox::Pointer<SparseDataType> dst);
   tbox::Pointer<SparseDataType> _createEmptySparseData();
   hier::Index _getRandomIndex();

   tbox::Pointer<SparseDataType> d_sparse_data;

   bool d_initialized;
   tbox::Dimension d_dim;
};

} // end namespace sam_test

#endif
