/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Test class for SparseData (implementation).
 *
 ************************************************************************/
#include "SparseDataTester.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputDatabase.h"

#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

#define NUM_INDICES 5
using namespace SAMRAI;

namespace sam_test {

SparseDataTester::SparseDataTester(
   const tbox::Dimension& dim):
   d_initialized(false),
   d_dim(dim)
{
}

SparseDataTester::~SparseDataTester() {
   d_sparse_data->clear();
}

bool
SparseDataTester::testConstruction()
{
   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));

   hier::IntVector vec(d_dim, 0);
   hier::IntVector ghosts(d_dim, 0);
   std::vector<std::string> dkeys;
   _getDblKeys(dkeys);
   std::vector<std::string> ikeys;
   _getIntKeys(ikeys);

   tbox::Pointer<SparseDataType> sparse(new SparseDataType(box, ghosts,
                                           dkeys, ikeys));
   d_sparse_data = sparse;

   d_sparse_data->printNames(tbox::plog);

   bool passed = true;
   if (!d_sparse_data->empty()) {
      passed = false;
      tbox::perr << "Sparse data should be empty and is not" << std::endl;
   } else {
      d_initialized = true;
   }

   SparseDataType::Iterator iter = d_sparse_data->registerIndex(hi);

   SparseDataType::AttributeIterator index_iter(d_sparse_data->begin(hi)),
   index_iterend(d_sparse_data->end(hi));

   if (index_iter != index_iterend) {
      passed = false;
      tbox::perr << "something is wrong with index2's list of attributes"
                 << std::endl;
   }

   d_sparse_data->clear();
   return passed;
}

bool
SparseDataTester::testCopy()
{

   //_fillObject(d_sparse_data);

   // ensure d_sparse_data is empty before we start
   d_sparse_data->clear();
   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();
   bool success = _testCopy(d_sparse_data, sample);
   // clean up d_sparse_data
   d_sparse_data->clear();
   sample->clear();
   return success;
}

bool
SparseDataTester::testCopy2()
{
   // ensure the tester's copy of d_sparse_data is empty before
   // we start
   d_sparse_data->clear();
   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();
   _fillObject(sample);

   bool success = _testCopy(sample, d_sparse_data);
   // clean up d_sparse_data
   d_sparse_data->clear();
   sample->clear();
   return success;
}

bool
SparseDataTester::testAdd()
{
   bool success = true;
   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();

   success = success && (sample->empty() ? true : false);
   if (success)
      _fillObject(sample);

   success = success && (!sample->empty() ? true : false);
   if (success)
      sample->clear();

   success = success && (sample->empty() ? true : false);
   return success;
}

bool
SparseDataTester::testRemove()
{
   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();
   _fillObject(sample);
   bool success = (!sample->empty() ? true : false);
   hier::Index idx = _getRandomIndex();

   int size = d_sparse_data->size();
   SparseDataType::Iterator iter = d_sparse_data->begin();
   while (iter != d_sparse_data->end()) {
      if (iter.getIndex() == idx) {
         sample->remove(iter);
      } else {
         ++iter;
      }
   }

   int newsize = sample->size();
   success = success && (newsize = size - 1) ? true : false;

   if (success) {
      sample->clear();
      success = success && (sample->empty() ? true : false);
   }

   return success;
}

bool
SparseDataTester::_testCopy(
   tbox::Pointer<SparseDataType> src,
   tbox::Pointer<SparseDataType> dst)
{
   src->copy(*dst);
   TBOX_ASSERT(src->size() == dst->size());
   pdat::SparseData<pdat::CellGeometry>::Iterator me(src);
   pdat::SparseData<pdat::CellGeometry>::Iterator me_end = src->end();
   pdat::SparseData<pdat::CellGeometry>::Iterator other(dst);

   bool success = true;
   for ( ; me != me_end && success != false; ++me, ++other) {
      if (me != other) {
         success = false;
      }
   }
   return success;
}

bool
SparseDataTester::testPackStream()
{
   bool success = true;

   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();
   _fillObject(sample);

   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));
   hier::BoxContainer blist(box);
   hier::Transformation trans(hier::IntVector::getZero(d_dim));

   pdat::CellOverlap overlap(blist, trans);

   int strsize = sample->getDataStreamSize(overlap);

   SparseDataType::Iterator iter(sample);

   tbox::MessageStream str(strsize, tbox::MessageStream::Write);
   sample->packStream(str, overlap);
   tbox::plog << "Printing sample1" << std::endl;
   sample->printAttributes(tbox::plog);

   tbox::Pointer<SparseDataType> sample2 = _createEmptySparseData();
   tbox::MessageStream upStr(strsize, tbox::MessageStream::Read);
   memcpy(upStr.getBufferStart(), str.getBufferStart(), strsize);

   sample2->unpackStream(upStr, overlap);

   SparseDataType::Iterator iter2(sample2);

   sample2->printNames(tbox::plog);
   tbox::plog << "Printing sample2" << std::endl;
   sample2->printAttributes(tbox::plog);
   for ( ; iter != sample->end() && iter2 != sample2->end(); ++iter, ++iter2) {
      tbox::plog << "iter1 node: " << std::endl;
      tbox::plog << iter;
      tbox::plog << "iter2 node: " << std::endl;
      tbox::plog << iter2;
      if (!iter.equals(iter2)) {
         success = false;
      }
      tbox::plog << std::endl;
   }

   const pdat::DoubleAttributeId did =
      sample2->getDblAttributeId("double_key_0");

   if (!sample2->isValid(did)) {
      success = false;
   }
   sample->clear();
   sample2->clear();
   return success;
}

bool
SparseDataTester::testDatabaseInterface()
{
   bool success = true;

   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();
   _fillObject(sample);
   tbox::Pointer<tbox::Database> input_db(new tbox::InputDatabase("input_db"));
   sample->putToDatabase(input_db);

   tbox::Pointer<SparseDataType> sample2 = _createEmptySparseData();
   sample2->getFromDatabase(input_db);

   SparseDataType::Iterator iter1(sample);
   SparseDataType::Iterator iter2(sample2);

   for ( ; iter1 != sample->end() && iter2 != sample2->end() && success;
         ++iter1, ++iter2) {
      if (!iter1.equals(iter2)) {
         success = false;
      }
   }

   sample->clear();
   sample2->clear();
   return success;
}

void
SparseDataTester::testTiming()
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->getTimer("SparseDataAddItem", true);

   hier::IntVector v(d_dim, 0);
   tbox::plog << "Begin Timing" << std::endl;
   tbox::Pointer<SparseDataType> sample = _createEmptySparseData();

   timer->start();
   SparseDataType::Iterator iter;
   for (int i = 0; i < 100; ++i) {
      for (int j = 0; j < 100; ++j) {
         v[0] = i;
         v[1] = j;
         hier::Index idx(v);
         double dvalues[DSIZE];
         _getDblValues(dvalues);
         int ivalues[ISIZE];
         _getIntValues(ivalues);

         iter = sample->registerIndex(idx);
         for (int m = i; m < j + 1; ++m) {
            iter.insert(dvalues, ivalues);
         }
      }
   }

   timer->stop();
   int numItems = sample->size();
   tbox::plog << numItems << std::endl;
   tbox::plog.precision(16);
   tbox::plog << "SparseData addItem insert time : "
              << timer->getTotalWallclockTime() << std::endl;
   tbox::plog << "End Timing" << std::endl;
   sample->clear();
}

hier::Index
SparseDataTester::_getRandomIndex()
{
   // return a random index in the range that we created
   // with _fillObject
   int val = (rand() % NUM_INDICES);
   hier::IntVector v(d_dim, 0);
   v[0] = val;
   v[1] = val;
   return hier::Index(v);
}

void
SparseDataTester::_fillObject(
   tbox::Pointer<SparseDataType> sparse_data)
{
   hier::IntVector v(d_dim, 0);
   double* dvalues = new double[DSIZE];
   _getDblValues(dvalues);
   int* ivalues = new int[ISIZE];
   _getIntValues(ivalues);

   SparseDataType::Iterator iter;
   for (int i = 1; i < NUM_INDICES; ++i) {
      v[0] = i;
      v[1] = i;
      hier::Index idx(v);
      iter = sparse_data->registerIndex(idx);
      for (int k = 0; k < i; ++k) {
         iter.insert(dvalues, ivalues);
      }
   }

   delete[] dvalues;
   delete[] ivalues;
}

tbox::Pointer<pdat::SparseData<pdat::CellGeometry> >
SparseDataTester::_createEmptySparseData()
{
   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));
   hier::IntVector ghosts(d_dim, 0);

   std::vector<std::string> dkeys;
   _getDblKeys(dkeys);
   std::vector<std::string> ikeys;
   _getIntKeys(ikeys);
   tbox::Pointer<SparseDataType> sample(
      new SparseDataType(box, ghosts, dkeys, ikeys));
   return sample;
}

void
SparseDataTester::_getDblKeys(std::vector<std::string>& keys)
{
   for (int i = 0; i < DSIZE; ++i) {
      std::stringstream key_name;
      key_name << "DOUBLE_KEY_" << i;
      keys.push_back(key_name.str());
   }
}

void
SparseDataTester::_getIntKeys(std::vector<std::string>& keys)
{
   for (int i = 0; i < ISIZE; ++i) {
      std::stringstream key_name;
      key_name << "INTEGER_KEY_" << i;
      keys.push_back(key_name.str());
   }
}

void
SparseDataTester::_getDblValues(double* values)
{
   for (int i = 0; i < DSIZE; ++i) {
      values[i] = (double)i;
   }
}

void
SparseDataTester::_getIntValues(int* values)
{
   for (int i = 0; i < ISIZE; ++i) {
      values[i] = i;
   }
}
} // end namespace sam_test
