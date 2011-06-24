/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence classes. 
 *
 ************************************************************************/

#ifndef included_xfer_RefineClasses_C
#define included_xfer_RefineClasses_C

#include <typeinfo>

#include "SAMRAI/xfer/RefineClasses.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/RefineClasses.I"
#endif

namespace SAMRAI {
namespace xfer {

int RefineClasses::s_default_refine_item_array_size = 20;

/*
 *************************************************************************
 *                                                                       *
 * Default constructor.                                                  *
 *                                                                       *
 *************************************************************************
 */

RefineClasses::RefineClasses():
   d_num_refine_items(0)
{
   d_equivalence_class_ids.resizeArray(0);
   d_refine_classes_data_items.resizeArray(s_default_refine_item_array_size);
}

/*
 *************************************************************************
 *									*
 * The destructor implicitly deletes the item storage associated with	*
 * the equivalence classes (and also the refine algorithm).		*
 *									*
 *************************************************************************
 */

RefineClasses::~RefineClasses()
{
   d_equivalence_class_ids.resizeArray(0);
   d_refine_classes_data_items.resizeArray(0);
}

/*
 *************************************************************************
 *                                                                       *
 * Return representative item for given equivalence class (first in list)*
 *                                                                       *
 *************************************************************************
 */

const RefineClasses::Data&
RefineClasses::getClassRepresentative(
   int equiv_class_id) const
{
   TBOX_ASSERT((equiv_class_id >= 0) &&
      (equiv_class_id < d_equivalence_class_ids.size()));
   return d_refine_classes_data_items[
             d_equivalence_class_ids[equiv_class_id].getFirstItem()];
}

/*
 *************************************************************************
 *                                                                       *
 * Return iterator for list of refine items for given equivalence class  *
 *                                                                       *
 *************************************************************************
 */

tbox::List<int>::Iterator
RefineClasses::getIterator(
   int equiv_class_id)
{
   TBOX_ASSERT((equiv_class_id >= 0) &&
      (equiv_class_id < d_equivalence_class_ids.size()));
   return tbox::List<int>::
          Iterator(d_equivalence_class_ids[equiv_class_id]);

}

/*
 *************************************************************************
 *									*
 * Insert a data item into the refine data list for the proper           *
 * equivalence class in sorted order by ascending operator priority.     *
 *									*
 *************************************************************************
 */

void RefineClasses::insertEquivalenceClassItem(
   RefineClasses::Data& data,
   tbox::Pointer<hier::PatchDescriptor> descriptor)
{

   if (!checkRefineItem(data, descriptor)) {
      tbox::perr << "Bad refine class data passed to "
                 << "RefineClasses::insertEquivalenceClassItem\n";
      printRefineItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = d_equivalence_class_ids.size();
         d_equivalence_class_ids.resizeArray(eq_index + 1);
      }

      data.d_class_id = eq_index;

      if (d_num_refine_items >= d_refine_classes_data_items.size()) {
         d_refine_classes_data_items.resizeArray(
            d_num_refine_items + s_default_refine_item_array_size);
      }

      d_refine_classes_data_items[d_num_refine_items] = data;

      d_equivalence_class_ids[eq_index].appendItem(d_num_refine_items);

      d_num_refine_items++;

   }

}

/*
 *************************************************************************
 *                                                                       *
 * Check for valid patch data ids, patch data types, and that scratch    *
 * data entry has at least as many ghost cells as destination data entry *
 * and stencil width of operator.  If so, return true; else return false.*
 * A descriptive error message is sent to TBOX_ERROR when a problem      *
 * appears.  If a null patch descriptor argument is passed, the          *
 * descriptor associated with the variable database Singleton object     *
 * will be used.                                                         *
 *                                                                       *
 *************************************************************************
 */

bool RefineClasses::checkRefineItem(
   const RefineClasses::Data& data_item,
   tbox::Pointer<hier::PatchDescriptor> descriptor) const
{

   bool item_good = true;

   tbox::Pointer<hier::PatchDescriptor> pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }

   const int dst_id = data_item.d_dst;
   const int src_id = data_item.d_src;
   const int scratch_id = data_item.d_scratch;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Source' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (scratch_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Scratch' patch data id invalid (< 0!)" << std::endl);
   }

   tbox::Pointer<hier::PatchDataFactory> dst_fact =
      pd->getPatchDataFactory(dst_id);
   tbox::Pointer<hier::PatchDataFactory> src_fact =
      pd->getPatchDataFactory(src_id);
   tbox::Pointer<hier::PatchDataFactory> scratch_fact =
      pd->getPatchDataFactory(scratch_id);

   if (item_good && !(src_fact->validCopyTo(scratch_fact))) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "It is not a valid operation to copy from `Source' patch data \n"
         << pd->mapIndexToName(src_id) << " to `Scratch' patch data "
         << pd->mapIndexToName(scratch_id) << std::endl);
   }

   if (item_good && !(scratch_fact->validCopyTo(dst_fact))) {
      item_good = false;
      pd->mapIndexToName(scratch_id);
      pd->mapIndexToName(dst_id);
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "It is not a valid operation to copy from `Scratch' patch data \n"
         << pd->mapIndexToName(scratch_id) << " to `Destination' patch data "
         << pd->mapIndexToName(dst_id) << std::endl);
   }

   const hier::IntVector& scratch_gcw = scratch_fact->getGhostCellWidth();

   if (item_good && (dst_fact->getGhostCellWidth() > scratch_gcw)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Destination' patch data " << pd->mapIndexToName(dst_id)
         << " has a larger ghost cell width than \n"
         << "`Scratch' patch data " << pd->mapIndexToName(scratch_id)
         << "\n`Destination' ghost width = "
         << dst_fact->getGhostCellWidth()
         << "\n`Scratch' ghost width = " << scratch_gcw << std::endl);
   }

   tbox::Pointer<hier::RefineOperator> refop = data_item.d_oprefine;
   if (item_good && !refop.isNull()) {
      if (refop->getStencilWidth() > scratch_gcw) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "Refine operator " << refop->getOperatorName()
            << "\nhas larger stencil width than ghost cell width"
            << "of `Scratch' patch data" << pd->mapIndexToName(scratch_id)
            << "\noperator stencil width = " << refop->getStencilWidth()
            << "\n`Scratch'  ghost width = " << scratch_gcw << std::endl);
      }
   }

   tbox::Pointer<VariableFillPattern> fill_pattern =
      data_item.d_var_fill_pattern;
   if (item_good && !fill_pattern.isNull()) {
      if (fill_pattern->getPatternName() != "BOX_GEOMETRY_FILL_PATTERN") {
         if (fill_pattern->getStencilWidth() > scratch_gcw) {
            item_good = false;
            TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
               << "VariableFillPattern " << fill_pattern->getPatternName()
               << "\nhas larger stencil width than ghost cell width"
               << "of `Scratch' patch data" << pd->mapIndexToName(
                  scratch_id)
               << "\noperator stencil width = "
               << fill_pattern->getStencilWidth()
               << "\n`Scratch'  ghost width = " << scratch_gcw << std::endl);
         }
      }
   }

   if (item_good && data_item.d_time_interpolate) {
      const int src_told_id = data_item.d_src_told;
      const int src_tnew_id = data_item.d_src_tnew;

      if (src_told_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source old' patch data id invalid (< 0!),\n"
            << "yet a request has made to time interpolate" << std::endl);
      }
      if (item_good && src_tnew_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source new' patch data id invalid (< 0!),\n"
            << "yet a request has made to time interpolate with them"
            << std::endl);
      }

      tbox::Pointer<hier::PatchDataFactory> src_told_fact =
         pd->getPatchDataFactory(src_told_id);
      tbox::Pointer<hier::PatchDataFactory> src_tnew_fact =
         pd->getPatchDataFactory(src_tnew_id);

      if (item_good && typeid(*src_told_fact) != typeid(*src_fact)) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " and `Source old' patch data "
            << pd->mapIndexToName(src_told_id)
            << " have different patch data types, yet a request has"
            << "\n been made to time interpolate with them" << std::endl);
      }

      if (item_good && typeid(*src_tnew_fact) != typeid(*src_fact)) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " and `Source new' patch data "
            << pd->mapIndexToName(src_tnew_id)
            << " have different patch data types, yet a request has"
            << "\n been made to time interpolate with them" << std::endl);
      }

   }

   return item_good;

}

/*
 *************************************************************************
 *                                                                       *
 * Compare refine data items in this refine classes object against       *
 * those in the argument refine classes object.  Return true if they     *
 * all match with regard to the patch data types, patch data ghost cell  *
 * widths, operator stencils, etc. that they refer to and return false   *
 * otherwise.  If a null patch descriptor argument is passed, the        *
 * descriptor associated with the variable database Singleton object     *
 * will be used.                                                         *
 *                                                                       *
 *************************************************************************
 */

bool RefineClasses::checkConsistency(
   tbox::Pointer<RefineClasses> test_classes,
   tbox::Pointer<hier::PatchDescriptor> descriptor) const
{

   bool items_match = true;

   tbox::Pointer<hier::PatchDescriptor> pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }

   if (d_equivalence_class_ids.size() !=
       test_classes->d_equivalence_class_ids.size()) {

      items_match = false;

   } else {

      int num_equiv_classes = d_equivalence_class_ids.size();
      int eq_index = 0;
      while (items_match && eq_index < num_equiv_classes) {

         if (d_equivalence_class_ids[eq_index].size() !=
             test_classes->
             d_equivalence_class_ids[eq_index].size()) {

            items_match = false;

         } else {

            tbox::List<int>::Iterator myli(d_equivalence_class_ids[eq_index]);
            tbox::List<int>::Iterator testli(test_classes->d_equivalence_class_ids[eq_index]);
            while (items_match && myli) {

               const RefineClasses::Data& myli_item =
                  d_refine_classes_data_items[myli()];
               const RefineClasses::Data& testli_item =
                  d_refine_classes_data_items[testli()];

               items_match = checkPatchDataItemConsistency(
                     myli_item.d_dst, testli_item.d_dst, pd);

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                        myli_item.d_src, testli_item.d_src, pd);
               }

               if (items_match && myli_item.d_time_interpolate) {
                  items_match = checkPatchDataItemConsistency(
                        myli_item.d_src_told,
                        testli_item.d_src_told,
                        pd);
               }

               if (items_match && myli_item.d_time_interpolate) {
                  items_match = checkPatchDataItemConsistency(
                        myli_item.d_src_tnew,
                        testli_item.d_src_tnew,
                        pd);
               }

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                        myli_item.d_scratch,
                        testli_item.d_scratch,
                        pd);
               }

               if (items_match) {
                  items_match = (myli_item.d_fine_bdry_reps_var ==
                                 testli_item.d_fine_bdry_reps_var);
               }

               if (items_match) {
                  items_match = (!myli_item.d_oprefine.isNull() ==
                                 !testli_item.d_oprefine.isNull());
                  if (items_match && !myli_item.d_oprefine.isNull()) {
                     items_match =
                        (myli_item.d_oprefine->getStencilWidth() ==
                         testli_item.d_oprefine->getStencilWidth());
                  }
               }

               if (items_match) {
                  items_match = (myli_item.d_time_interpolate ==
                                 testli_item.d_time_interpolate);
                  if (items_match && myli_item.d_time_interpolate) {
                     items_match = (typeid(*(myli_item.d_optime)) ==
                                    typeid(*(testli_item.d_optime)));
                  }
               }

               if (items_match) {
                  items_match = (!myli_item.d_var_fill_pattern.isNull() ==
                                 !testli_item.d_var_fill_pattern.isNull());
                  if (items_match && !myli_item.d_var_fill_pattern.isNull()) {
                     items_match = (typeid(*(myli_item.d_var_fill_pattern)) ==
                                    typeid(*(testli_item.d_var_fill_pattern)));
                  }
               }

               myli++;
               testli++;

            } // while items in equivalence class match

         } // if number of items in equivalence class match

         eq_index++;

      } // while equivalence classes match

   } // else number of equivalence classes do not match

   return items_match;

}

/*
 *************************************************************************
 *                                                                       *
 * Private member function to determine whether two patch data items     *
 * are consistent.                                                       *
 *                                                                       *
 *************************************************************************
 */

bool RefineClasses::checkPatchDataItemConsistency(
   int item_id1,
   int item_id2,
   tbox::Pointer<hier::PatchDescriptor> pd) const
{

   bool items_match = ((item_id1 >= 0) && (item_id2 >= 0));

   if (items_match) {

      tbox::Pointer<hier::PatchDataFactory> pdf1 =
         pd->getPatchDataFactory(item_id1);
      tbox::Pointer<hier::PatchDataFactory> pdf2 =
         pd->getPatchDataFactory(item_id2);

      items_match = (typeid(*pdf1) == typeid(*pdf2));

      if (items_match) {
         items_match = (pdf1->getGhostCellWidth() ==
                        pdf2->getGhostCellWidth());
      }

   }

   return items_match;

}

/*
 *************************************************************************
 *                                                                       *
 * Private member function to determine equivalence class for            *
 * given data item.  Return value of -1 indicates no match found; else   *
 * return value is index of match.                                       *
 *                                                                       *
 *************************************************************************
 */

int RefineClasses::getEquivalenceClassIndex(
   const RefineClasses::Data& data,
   tbox::Pointer<hier::PatchDescriptor> descriptor) const
{

   int eq_index = -1;

   tbox::Pointer<hier::PatchDescriptor> pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }

   int dst_id = data.d_dst;
   int src_id = data.d_src;

   tbox::Pointer<hier::PatchDataFactory> dst_pdf =
      pd->getPatchDataFactory(dst_id);
   tbox::Pointer<hier::PatchDataFactory> src_pdf =
      pd->getPatchDataFactory(src_id);

   hier::IntVector dst_ghosts = dst_pdf->getGhostCellWidth();
   hier::IntVector src_ghosts = src_pdf->getGhostCellWidth();

   int num_equiv_classes = d_equivalence_class_ids.size();
   bool equiv_found = false;
   int nl = 0;
   while (!equiv_found && nl < num_equiv_classes) {

      bool dst_equiv = false;
      bool src_equiv = false;

      const RefineClasses::Data& class_rep = getClassRepresentative(nl);

      int rep_dst_id = class_rep.d_dst;
      tbox::Pointer<hier::PatchDataFactory> rep_dst_pdf =
         pd->getPatchDataFactory(rep_dst_id);
      hier::IntVector rep_dst_ghosts =
         rep_dst_pdf->getGhostCellWidth();

      /*
       * Check if destinations are equivalent
       */
      if ((dst_ghosts == rep_dst_ghosts) &&
          (typeid(*dst_pdf) == typeid(*rep_dst_pdf))) {
         dst_equiv = true;
      }

      /*
       * If src_id and dst_id are the same, there is nothing more to check.
       * Otherwise, if destinations were equivalent, check if sources
       * are equivalent.
       */
      if (dst_id == src_id) {
         src_equiv = dst_equiv;
      } else if (dst_equiv) {
         int rep_src_id = class_rep.d_src;
         tbox::Pointer<hier::PatchDataFactory> rep_src_pdf =
            pd->getPatchDataFactory(rep_src_id);
         hier::IntVector rep_src_ghosts =
            rep_src_pdf->getGhostCellWidth();
         if ((src_ghosts == rep_src_ghosts) &&
             (typeid(*src_pdf) == typeid(*rep_src_pdf))) {
            src_equiv = true;
         }
      }

      /*
       * Check if fill patterns are the same.
       */
      bool fill_patterns_same = false;
      if (dst_equiv && src_equiv) {
         fill_patterns_same =
            (data.d_var_fill_pattern->getPatternName() ==
             class_rep.d_var_fill_pattern->getPatternName());
      }

      /*
       * If destinations and sources are both equivalent, exit loop
       * and set return value to identify current equivalence class id.
       */
      if (dst_equiv && src_equiv && fill_patterns_same) {
         eq_index = nl;
         equiv_found = true;
      }

      nl++;

   }

   return eq_index;

}

/*
 *************************************************************************
 *                                                                       *
 * Increase the data items array to the specified size.                  *
 *                                                                       *
 *************************************************************************
 */

void RefineClasses::increaseRefineItemArraySize(
   const int size)
{
   if (size > d_refine_classes_data_items.size()) {
      d_refine_classes_data_items.resizeArray(size);
   }
}

/*
 *************************************************************************
 *									*
 * Print the data in the refine item lists to the specified stream.      *
 *									*
 *************************************************************************
 */

void RefineClasses::printClassData(
   std::ostream& stream) const
{
   stream << "RefineClasses::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < d_equivalence_class_ids.size(); i++) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      for (tbox::List<int>::Iterator
           li(d_equivalence_class_ids[i]); li; li++) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printRefineItem(stream, d_refine_classes_data_items[li()]);

         j++;
      }
      stream << std::endl;
   }

}

void RefineClasses::printRefineItem(
   std::ostream& stream,
   const RefineClasses::Data& data) const
{
   stream << "\n";
   stream << "desination component:   "
          << data.d_dst << std::endl;
   stream << "source component:       "
          << data.d_src << std::endl;
   stream << "scratch component:      "
          << data.d_scratch << std::endl;
   stream << "fine boundary represents variable:      "
          << data.d_fine_bdry_reps_var << std::endl;
   stream << "tag:      "
          << data.d_tag << std::endl;

   if (data.d_oprefine.isNull()) {
      stream << "NULL refining operator" << std::endl;
   } else {
      stream << "refine operator name:          "
             << typeid(*data.d_oprefine).name()
             << std::endl;
      stream << "operator priority:      "
             << data.d_oprefine->getOperatorPriority()
             << std::endl;
      stream << "operator stencil width: "
             << data.d_oprefine->getStencilWidth()
             << std::endl;
   }
   if (!data.d_time_interpolate) {
      stream << "time interpolate is false" << std::endl;
   } else {
      stream << "old source component:   "
             << data.d_src_told << std::endl;
      stream << "new source component:   "
             << data.d_src_tnew << std::endl;
      stream << "time interpolation operator name:          "
             << typeid(*data.d_optime).name()
             << std::endl;
   }
   if (data.d_var_fill_pattern.isNull()) {
      stream << "var fill pattern is null" << std::endl;
   } else {
      stream << "var fill pattern name:          "
             << typeid(*data.d_var_fill_pattern).name()
             << std::endl;
   }
   stream << std::endl;
}

}
}
#endif
