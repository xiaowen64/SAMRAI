/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes. 
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenClasses_C
#define included_xfer_CoarsenClasses_C

#include "SAMRAI/xfer/CoarsenClasses.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/CoarsenClasses.I"
#endif

#include <typeinfo>

namespace SAMRAI {
namespace xfer {

int CoarsenClasses::s_default_coarsen_item_array_size = 20;

/*
 *************************************************************************
 *                                                                       *
 * Constructor sets boolean for filling coarse data and creates new      *
 * array of equivalence classes.                                         *
 *                                                                       *
 *************************************************************************
 */

CoarsenClasses::CoarsenClasses(
   bool fill_coarse_data):
   d_fill_coarse_data(fill_coarse_data),
   d_num_coarsen_items(0)
{
   d_coarsen_classes_data_items.resizeArray(
      s_default_coarsen_item_array_size,
      Data(tbox::Dimension::getInvalidDimension()));
}

/*
 *************************************************************************
 *									*
 * The destructor implicitly deletes the item storage associated with	*
 * the equivalence classes (and also the coarsen algorithm).		*
 *									*
 *************************************************************************
 */

CoarsenClasses::~CoarsenClasses()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Return representative item for given equivalence class (first in list)*
 *                                                                       *
 *************************************************************************
 */

const CoarsenClasses::Data&
CoarsenClasses::getClassRepresentative(
   int equiv_class_id) const
{
   TBOX_ASSERT((equiv_class_id >= 0) &&
      (equiv_class_id < (int)d_equivalence_class_ids.size()));

   return d_coarsen_classes_data_items[
             d_equivalence_class_ids[equiv_class_id].getFirstItem()];
}

/*
 *************************************************************************
 *                                                                       *
 * Return iterator for list of coarsen items for given equivalence class  *
 *                                                                       *
 *************************************************************************
 */

tbox::List<int>::Iterator
CoarsenClasses::getIterator(
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
 * Insert a data item into the coarsen data list for the proper          *
 * equivalence class in sorted order by ascending operator priority.     *
 *									*
 *************************************************************************
 */

void CoarsenClasses::insertEquivalenceClassItem(
   CoarsenClasses::Data& data,
   tbox::Pointer<hier::PatchDescriptor> descriptor)
{

   if (!checkCoarsenItem(data, descriptor)) {
      tbox::perr << "Bad coarsen class data passed to "
                 << "CoarsenClasses::insertEquivalenceClassItem\n";
      printCoarsenItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = d_equivalence_class_ids.size();
         d_equivalence_class_ids.resizeArray(eq_index + 1);
      }

      data.d_class_id = eq_index;

      if (d_num_coarsen_items >= d_coarsen_classes_data_items.size()) {
         d_coarsen_classes_data_items.resizeArray(
            d_num_coarsen_items + s_default_coarsen_item_array_size,
            Data(data.d_gcw_to_coarsen.getDim()));
      }

      d_coarsen_classes_data_items[d_num_coarsen_items] = data;

      d_equivalence_class_ids[eq_index].appendItem(d_num_coarsen_items);

      d_num_coarsen_items++;
   }

}

/*
 *************************************************************************
 *                                                                       *
 * Check for valid patch data ids, patch data types, and that source and *
 * destination data entries have sufficient ghost cells to satisfy the   *
 * coarsen operator and necessary copy operations.  If so, return true;  *
 * else return false.  A descriptive error message is sent to TBOX_ERROR *
 * when a problem appears.  If a null patch descriptor argument is       *
 * passed, the descriptor associated with the variable database          *
 * Singleton object will be used.                                        *
 *                                                                       *
 *************************************************************************
 */

bool CoarsenClasses::checkCoarsenItem(
   const CoarsenClasses::Data& data_item,
   tbox::Pointer<hier::PatchDescriptor> descriptor) const
{

   bool item_good = true;

   tbox::Pointer<hier::PatchDescriptor> pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }

   const int dst_id = data_item.d_dst;
   const int src_id = data_item.d_src;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "`Source' patch data id invalid (< 0!)" << std::endl);
   }

   tbox::Pointer<hier::PatchDataFactory> dfact =
      pd->getPatchDataFactory(dst_id);
   tbox::Pointer<hier::PatchDataFactory> sfact =
      pd->getPatchDataFactory(src_id);

   if (item_good && !(sfact->validCopyTo(dfact))) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "It is not a valid operation to copy from `Source' patch data \n"
         << pd->mapIndexToName(src_id) << " to `Destination' patch data "
         << pd->mapIndexToName(dst_id) << std::endl);
   }

   tbox::Pointer<hier::CoarsenOperator> coarsop = data_item.d_opcoarsen;
   if (item_good && !coarsop.isNull()) {
      if (coarsop->getStencilWidth() > sfact->getGhostCellWidth()) {
         item_good = false;
         TBOX_ERROR("Bad data given to CoarsenClasses...\n"
            << "Coarsen operator " << coarsop->getOperatorName()
            << "\nhas larger stencil width than ghost cell width"
            << "of `Source' patch data" << pd->mapIndexToName(src_id)
            << "\noperator stencil width = " << coarsop->getStencilWidth()
            << "\n`Source'  ghost width = "
            << sfact->getGhostCellWidth()
            << std::endl);
      }
   }

   return item_good;

}

/*
 *************************************************************************
 *                                                                       *
 * Compare coarsen data items in this coarsen classes object against     *
 * those in the argument coarsen classes object.  Return true if they    *
 * all match with regard to the patch data types, patch data ghost cell  *
 * widths, operator stencils, etc. that they refer to and return false   *
 * otherwise.  If a null patch descriptor argument is passed, the        *
 * descriptor associated with the variable database Singleton object     *
 * will be used.                                                         *
 *                                                                       *
 *************************************************************************
 */

bool CoarsenClasses::checkConsistency(
   tbox::Pointer<CoarsenClasses> test_classes,
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
            tbox::List<int>::Iterator testli(d_equivalence_class_ids[eq_index]);
            while (items_match && myli) {

               const CoarsenClasses::Data& myli_item =
                  d_coarsen_classes_data_items[myli()];
               const CoarsenClasses::Data& testli_item =
                  d_coarsen_classes_data_items[testli()];

               items_match = checkPatchDataItemConsistency(
                     myli_item.d_dst, testli_item.d_dst, pd);

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                        myli_item.d_src, testli_item.d_src, pd);
               }

               if (items_match) {
                  items_match = (myli_item.d_fine_bdry_reps_var ==
                                 testli_item.d_fine_bdry_reps_var);
               }

               if (items_match) {
                  items_match = (myli_item.d_gcw_to_coarsen ==
                                 testli_item.d_gcw_to_coarsen);
               }

               if (items_match) {
                  items_match = (!myli_item.d_opcoarsen.isNull() ==
                                 !testli_item.d_opcoarsen.isNull());
                  if (items_match && !myli_item.d_opcoarsen.isNull()) {
                     items_match =
                        (myli_item.d_opcoarsen->getStencilWidth() ==
                         testli_item.d_opcoarsen->getStencilWidth());
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

bool CoarsenClasses::checkPatchDataItemConsistency(
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

int CoarsenClasses::getEquivalenceClassIndex(
   const CoarsenClasses::Data& data,
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

      const CoarsenClasses::Data& class_rep = getClassRepresentative(nl);

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
 * Return the number of items in the specified equivalence class.        *
 *                                                                       *
 *************************************************************************
 */

int CoarsenClasses::getNumberOfItemsInEquivalenceClass(
   int equiv_class_id) const
{
   return d_equivalence_class_ids[equiv_class_id].size();
}

/*
 *************************************************************************
 *                                                                       *
 * Increase the data items array to the specified size.                  *
 *                                                                       *
 *************************************************************************
 */

void CoarsenClasses::increaseCoarsenItemArraySize(
   const int size,
   const tbox::Dimension& dim)
{
   if (size > d_coarsen_classes_data_items.size()) {
      d_coarsen_classes_data_items.resizeArray(size, Data(dim));
   }
}

/*
 *************************************************************************
 *									*
 * Print the data in the coarsen item lists to the specified stream.     *
 *									*
 *************************************************************************
 */

void CoarsenClasses::printClassData(
   std::ostream& stream) const
{
   stream << "CoarsenClasses::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < (int)d_equivalence_class_ids.size(); i++) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      for (tbox::List<int>::Iterator
           li(d_equivalence_class_ids[i]); li; li++) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printCoarsenItem(stream, d_coarsen_classes_data_items[li()]);

         j++;
      }
      stream << std::endl;
   }

}

void CoarsenClasses::printCoarsenItem(
   std::ostream& stream,
   const CoarsenClasses::Data& data) const
{
   stream << "\n";
   stream << "desination component:   "
          << data.d_dst << std::endl;
   stream << "source component:       "
          << data.d_src << std::endl;
   stream << "fine boundary represents variable:       "
          << data.d_fine_bdry_reps_var << std::endl;
   stream << "gcw to coarsen:       "
          << data.d_gcw_to_coarsen << std::endl;
   stream << "tag:       "
          << data.d_tag << std::endl;

   if (data.d_opcoarsen.isNull()) {
      stream << "NULL coarsening operator" << std::endl;
   } else {
      stream << "coarsen operator name:          "
             << typeid(*data.d_opcoarsen).name()
             << std::endl;
      stream << "operator priority:      "
             << data.d_opcoarsen->getOperatorPriority()
             << std::endl;
      stream << "operator stencil width: "
             << data.d_opcoarsen->getStencilWidth()
             << std::endl;
   }
   stream << std::endl;
}

}
}
#endif
