/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test cell-centered patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
using namespace std;

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyCellDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include "boost/shared_ptr.hpp"

using namespace SAMRAI;

/* Helper function declarations */
bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy);

#define NVARS 4

int main(
   int argc,
   char* argv[]) {

   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI::MAX_DIM_VAL);
   const tbox::Dimension dim(d);

   const std::string log_fn = std::string("cell_hieops.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      int ln, iv;

      /*
       * Make a simple 2-level hierarchy.
       */
      double lo[SAMRAI::MAX_DIM_VAL];
      double hi[SAMRAI::MAX_DIM_VAL];

      hier::Index clo0(dim);
      hier::Index chi0(dim);
      hier::Index clo1(dim);
      hier::Index chi1(dim);
      hier::Index flo0(dim);
      hier::Index fhi0(dim);
      hier::Index flo1(dim);
      hier::Index fhi1(dim);

      for (int i = 0; i < dim.getValue(); i++) {
         lo[i] = 0.0;
         clo0(i) = 0;
         flo0(i) = 4;
         fhi0(i) = 7;
         if (i == 1) {
            hi[i] = 0.5;
            chi0(i) = 2;
            clo1(i) = 3;
            chi1(i) = 4;
         } else {
            hi[i] = 1.0;
            chi0(i) = 9;
            clo1(i) = 0;
            chi1(i) = 9;
         }
         if (i == 0) {
            flo1(i) = 8;
            fhi1(i) = 13;
         } else {
            flo1(i) = flo0(i);
            fhi1(i) = fhi0(i);
         }
      }

      hier::Box coarse0(clo0, chi0, hier::BlockId(0));
      hier::Box coarse1(clo1, chi1, hier::BlockId(0));
      hier::Box fine0(flo0, fhi0, hier::BlockId(0));
      hier::Box fine1(flo1, fhi1, hier::BlockId(0));
      hier::IntVector ratio(dim, 2);

      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_domain;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_domain.pushBack(fine0);
      fine_domain.pushBack(fine1);

      boost::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      const int nproc = mpi.getSize();

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_domain.size();

      hier::BoxLevel layer0(hier::IntVector(dim, 1), geometry);
      hier::BoxLevel layer1(ratio, geometry);

      hier::BoxContainer::iterator coarse_itr(coarse_domain);
      for (int ib = 0; ib < n_coarse_boxes; ib++, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0.getMPI().getRank()) {
               layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib),
                  layer0.getMPI().getRank()));
            }
         } else {
            layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_itr(fine_domain);
      for (int ib = 0; ib < n_fine_boxes; ib++, ++fine_itr) {
         if (nproc > 1) {
            if (ib == layer1.getMPI().getRank()) {
               layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib),
                  layer1.getMPI().getRank()));
            }
         } else {
            layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      /*
       * Create some variables, a context, and register them with
       * the variable database.
       */
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      boost::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      boost::shared_ptr<pdat::CellVariable<double> > cvar[NVARS];
      int cvindx[NVARS];
      cvar[0].reset(new pdat::CellVariable<double>(dim, "cvar0", 1));
      cvindx[0] = variable_db->registerVariableAndContext(
            cvar[0], dummy, no_ghosts);
      cvar[1].reset(new pdat::CellVariable<double>(dim, "cvar1", 1));
      cvindx[1] = variable_db->registerVariableAndContext(
            cvar[1], dummy, no_ghosts);
      cvar[2].reset(new pdat::CellVariable<double>(dim, "cvar2", 1));
      cvindx[2] = variable_db->registerVariableAndContext(
            cvar[2], dummy, no_ghosts);
      cvar[3].reset(new pdat::CellVariable<double>(dim, "cvar3", 1));
      cvindx[3] = variable_db->registerVariableAndContext(
            cvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::CellVariable<double> > cwgt(
         new pdat::CellVariable<double>(dim, "cwgt", 1));
      int cwgt_id = variable_db->registerVariableAndContext(
            cwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(cvindx[iv]);
         }
      }

      /*
       * Create instances of hierarchy operations to apply certain
       * mathematical operators.  e.g. scale(), axpy(), min(), etc.
       */
      int coarsest = 0;
      int finest = 1;
      boost::shared_ptr<math::HierarchyDataOpsReal<double> > cell_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));
      TBOX_ASSERT(cell_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > cwgt_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));

      boost::shared_ptr<hier::Patch> patch;

      // Initialize control volume data for cell-centered components
      hier::Box coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
               patch->getPatchGeometry(),
               boost::detail::dynamic_cast_tag());
            const double* dx = pgeom->getDx();
            double cell_vol = dx[0];
            for (int i = 1; i < dim.getValue(); i++) {
               cell_vol *= dx[i];
            }
            boost::shared_ptr<pdat::CellData<double> > cvdata(
               patch->getPatchData(cwgt_id),
               boost::detail::dynamic_cast_tag());
            cvdata->fillAll(cell_vol);
            if (ln == 0) cvdata->fillAll(0.0, (coarse_fine * patch->getBox()));
         }
      }

      /*
       * Apply various operations to the hierarchy data. Test the
       * result to assure its accuracy.
       */
      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level
      bool vol_test_passed = true;
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            boost::shared_ptr<pdat::CellData<double> > cvdata(
               patch->getPatchData(cwgt_id),
               boost::detail::dynamic_cast_tag());

            pdat::CellIterator cend(cvdata->getBox(), false);
            for (pdat::CellIterator c(cvdata->getBox(), true);
                 c != cend && vol_test_passed; ++c) {
               pdat::CellIndex cell_index = *c;

               if (ln == 0) {
                  if ((coarse_fine * patch->getBox()).contains(cell_index)) {
                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), 0.0)) {
                        vol_test_passed = false;
                     }
                  } else {
                     double compare;
                     if ((dim == tbox::Dimension(2))) {
                        compare = 0.01;
                     }
                     else {
                        compare = 0.001;
                     }

                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), compare)) {

                        vol_test_passed = false;
                     }
                  }
               }

               if (ln == 1) {
                  double compare;
                  if ((dim == tbox::Dimension(2))) {
                     compare = 0.0025;
                  }
                  else {
                     compare = 0.000125;
                  }

                  if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                            cell_index), compare)) {
                     vol_test_passed = false;
                  }
               }
            }
         }
      }
      if (!vol_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1a: Check control volume data set properly"
         << endl;
         cwgt_ops->printData(cwgt_id, tbox::plog);
      }

      // Test #1b: HierarchyCellDataOpsReal2::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1b: HierarchyCellDataOpsReal2::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << endl;
      }

      // Test #2: HierarchyCellDataOpsReal2::numberOfEntries()
      // Expected: num_data_points = 90
      int num_data_points = cell_ops->numberOfEntries(cvindx[0]);

      {
         int compare;
         if ((dim == tbox::Dimension(2))) {
            compare = 90;
         }
         else {
            compare = 660;
         }

         if (num_data_points != compare) {
            num_failures++;
            tbox::perr
            << "FAILED: - Test #2: math::HierarchyCellDataOpsReal::numberOfEntries()\n"
            << "Expected value = " << compare << ", Computed value = "
            << num_data_points << std::endl;
         }
      }

      // Test #3a: HierarchyCellDataOpsReal2::setToScalar()
      // Expected: v0 = 2.0
      double val0 = 2.0;
      cell_ops->setToScalar(cvindx[0], val0);
      if (!doubleDataSameAsValue(cvindx[0], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3a: HierarchyCellDataOpsReal2::setToScalar()\n"
         << "Expected: v0 = " << val0 << endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #3b: HierarchyCellDataOpsReal2::setToScalar()
      // Expected: v1 = (4.0)
      cell_ops->setToScalar(cvindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(cvindx[1], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3b: HierarchyCellDataOpsReal2::setToScalar()\n"
         << "Expected: v1 = " << val1 << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #4: HierarchyCellDataOpsReal2::copyData()
      // Expected: v2 = v1 = (4.0)
      cell_ops->copyData(cvindx[2], cvindx[1]);
      if (!doubleDataSameAsValue(cvindx[2], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4: HierarchyCellDataOpsReal2::copyData()\n"
         << "Expected: v2 = " << val1 << endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #5: HierarchyCellDataOpsReal2::swapData()
      // Expected: v0 = (4.0), v1 = (2.0)
      cell_ops->swapData(cvindx[0], cvindx[1]);
      if (!doubleDataSameAsValue(cvindx[0], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5a: HierarchyCellDataOpsReal2::swapData()\n"
         << "Expected: v0 = " << val1 << endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(cvindx[1], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5b: HierarchyCellDataOpsReal2::swapData()\n"
         << "Expected: v1 = " << val0 << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #6: HierarchyCellDataOpsReal2::scale()
      // Expected: v2 = 0.25 * v2 = (1.0)
      cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_scale, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #6: HierarchyCellDataOpsReal2::scale()\n"
         << "Expected: v2 = " << val_scale << endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #7: HierarchyCellDataOpsReal2::add()
      // Expected: v3 = v0 + v1 = (6.0)
      cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(cvindx[3], val_add, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #7: HierarchyCellDataOpsReal2::add()\n"
         << "Expected: v3 = " << val_add << endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      cell_ops->setToScalar(cvindx[0], 0.0);

      // Test #8: HierarchyCellDataOpsReal2::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_sub, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #8: HierarchyCellDataOpsReal2::subtract()\n"
         << "Expected: v1 = " << val_sub << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9a: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v1 = v1 + (0.0) = (6.0)
      cell_ops->addScalar(cvindx[1], cvindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9a: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9b: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v2 = v2 + (0.0) = (1.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 0.0);
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9b: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #9c: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v2 = v2 + (3.0) = (4.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9c: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Reset v3:  v3 = (0.5)
      cell_ops->setToScalar(cvindx[3], 0.5);

      // Test #10: HierarchyCellDataOpsReal2::multiply()
      // Expected: v1 = v3 * v1 = (3.0)
      cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(cvindx[1], val_mult, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #10: HierarchyCellDataOpsReal2::multiply()\n"
         << "Expected: v1 = " << val_mult << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #11: HierarchyCellDataOpsReal2::divide()
      // Expected: v0 = v2 / v1 = 1.3333333333
      cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
      double val_div = 1.33333333333;
      if (!doubleDataSameAsValue(cvindx[0], val_div, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #11: HierarchyCellDataOpsReal2::divide()\n"
         << "Expected: v0 = " << val_div << endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #12: HierarchyCellDataOpsReal2::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.333333333)
      cell_ops->reciprocal(cvindx[1], cvindx[1]);
      double val_rec = 0.33333333333;
      if (!doubleDataSameAsValue(cvindx[1], val_rec, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #12: HierarchyCellDataOpsReal2::reciprocal()\n"
         << "Expected: v1 = " << val_rec << endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #13: HierarchyCellDataOpsReal2::abs()
      // Expected:  v3 = abs(v2) = 4.0
      cell_ops->abs(cvindx[3], cvindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(cvindx[3], val_abs, hierarchy)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #13: HierarchyCellDataOpsReal2::abs()\n"
                    << "Expected: v3 = " << val_abs << endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #14: Place some bogus values on coarse level
      boost::shared_ptr<pdat::CellData<double> > cdata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
         cdata = boost::dynamic_pointer_cast<pdat::CellData<double>,
                                             hier::PatchData>(patch->getPatchData(cvindx[2]));
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::CellIndex(index0), 0) = 100.0;
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::CellIndex(index1), 0) = -1000.0;
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level->begin());
           ipp != level->end(); ++ipp) {
         patch = *ipp;
         cdata = boost::dynamic_pointer_cast<pdat::CellData<double>,
                                             hier::PatchData>(patch->getPatchData(cvindx[2]));
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;

         pdat::CellIterator cend(cdata->getBox(), false);
         for (pdat::CellIterator c(cdata->getBox(), true);
              c != cend && bogus_value_test_passed; ++c) {
            pdat::CellIndex cell_index = *c;

            if (cell_index == pdat::CellIndex(index0)) {
               if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index),
                      100.0)) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (cell_index == pdat::CellIndex(index1)) {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), -1000.0)) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), 4.0)) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #14:  Place some bogus values on coarse level"
         << endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #15: HierarchyCellDataOpsReal2::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1452
      double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
      double compare;
      if ((dim == tbox::Dimension(2))) {
         compare = 1452;
      }
      else {
         compare = 3732;
      }
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, compare)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #15: HierarchyCellDataOpsReal2::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = " << compare << ", Computed value = "
         << setprecision(12) << bogus_l1_norm << endl;
      }

      // Test #16: HierarchyCellDataOpsReal2::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.0
      double correct_l1_norm = cell_ops->L1Norm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #16: HierarchyCellDataOpsReal2::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.0, Computed value = "
         << correct_l1_norm << endl;
      }

      // Test #17: HierarchyCellDataOpsReal2::L2Norm()
      // Expected:  l2_norm = 2.82842712475
      double l2_norm = cell_ops->L2Norm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 2.82842712475)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #17: HierarchyCellDataOpsReal2::L2Norm()\n"
         << "Expected value = 2.82842712475, Computed value = "
         << l2_norm << endl;
      }

      // Test #18: HierarchyCellDataOpsReal2::L2Norm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.0
      double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #18: HierarchyCellDataOpsReal2::L2Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << endl;
      }

      // Test #19: HierarchyCellDataOpsReal2::L2Norm() - w/control weight
      // Expected:  max_norm = 4.0
      double max_norm = cell_ops->maxNorm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19: HierarchyCellDataOpsReal2::L2Norm()"
         << " - w/control weight\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << endl;
      }

      // Reset data and test sums, axpy's
      cell_ops->setToScalar(cvindx[0], 1.00);
      cell_ops->setToScalar(cvindx[1], 2.5);
      cell_ops->setToScalar(cvindx[2], 7.0);

      // Test #20: HierarchyCellDataOpsReal2::linearSum()
      // Expected:  v3 = 5.0
      cell_ops->linearSum(cvindx[3], 2.0, cvindx[1], 0.00, cvindx[0]);
      double val_linearSum = 5.0;
      if (!doubleDataSameAsValue(cvindx[3], val_linearSum, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #20: HierarchyCellDataOpsReal2::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #21: HierarchyCellDataOpsReal2::axmy()
      // Expected:  v3 = 6.5
      cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(cvindx[3], val_axmy, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21: HierarchyCellDataOpsReal2::axmy()\n"
         << "Expected: v3 = " << val_axmy << endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #22a: HierarchyCellDataOpsReal2::dot() - (ind2) * (ind1)
      // Expected:  cdot = 8.75
      double cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22a: HierarchyCellDataOpsReal2::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
      }

      // Test #22b: HierarchyCellDataOpsReal2::dot() - (ind1) * (ind2)
      // Expected:  cdot = 8.75
      cdot = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22b: HierarchyCellDataOpsReal2::dot() - (ind1) * (ind2)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(cvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; iv++) {
         cvar[iv].reset();
      }
      cwgt.reset();

      geometry.reset();
      hierarchy.reset();
      cell_ops.reset();
      cwgt_ops.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  cell hierops" << std::endl;
      }

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   bool test_passed = true;

   int ln;
   boost::shared_ptr<hier::Patch> patch;
   for (ln = 0; ln < 2; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
         boost::shared_ptr<pdat::CellData<double> > cvdata(
            patch->getPatchData(desc_id),
            boost::detail::dynamic_cast_tag());

         pdat::CellIterator cend(cvdata->getBox(), false);
         for (pdat::CellIterator c(cvdata->getBox(), true);
              c != cend && test_passed; ++c) {
            pdat::CellIndex cell_index = *c;
            if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
