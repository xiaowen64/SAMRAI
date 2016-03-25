//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/standard/RefineClasses.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description:	Simple structure for managing refinement data in equivalence classes.
//

#ifndef included_xfer_RefineClasses_C
#define included_xfer_RefineClasses_C

#include <typeinfo>

#include "RefineClasses.h"

#include "IntVector.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "VariableDatabase.h"
#include "tbox/Utilities.h"


namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor.                                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM>  RefineClasses<DIM>::RefineClasses()
{
   d_refine_equivalence_classes.resizeArray(0);
}

/*
*************************************************************************
*									*
* The destructor implicitly deletes the item storage associated with	*
* the equivalence classes (and also the refine algorithm).		*
*									*
*************************************************************************
*/

template<int DIM>  RefineClasses<DIM>::~RefineClasses()
{
   d_refine_equivalence_classes.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Return number of equivalence classes.                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> int RefineClasses<DIM>::getNumberEquivalenceClasses() const
{
   return(d_refine_equivalence_classes.getSize());
}

/*
*************************************************************************
*                                                                       *
* Return representative item for given equivalence class (first in list)*
*                                                                       *
*************************************************************************
*/

template<int DIM> const typename RefineClasses<DIM>::Data& 
   RefineClasses<DIM>::getClassRepresentative(int equiv_class_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (equiv_class_id >= 0) && 
           (equiv_class_id < d_refine_equivalence_classes.getSize()) );
#endif 
   return(d_refine_equivalence_classes[equiv_class_id].getFirstItem());
}

/*
*************************************************************************
*                                                                       *
* Return iterator for list of refine items for given equivalence class  *
*                                                                       *
*************************************************************************
*/

template<int DIM> typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator 
RefineClasses<DIM>::getIterator(int equiv_class_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (equiv_class_id >= 0) && 
           (equiv_class_id < d_refine_equivalence_classes.getSize()) );
#endif 
   return(typename tbox::List<typename RefineClasses<DIM>::Data>::
          Iterator(d_refine_equivalence_classes[equiv_class_id]));
}

/*
*************************************************************************
*									*
* Insert a data item into the refine data list for the proper           *
* equivalence class in sorted order by ascending operator priority.     *
*									*
*************************************************************************
*/

template<int DIM> void RefineClasses<DIM>::insertEquivalenceClassItem(
   const typename RefineClasses<DIM>::Data& data,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor)
{

   if (!checkRefineItem(data, descriptor)) {
      tbox::perr << "Bad refine class data passed to "
           << "RefineClasses<DIM>::insertEquivalenceClassItem\n";
      printRefineItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = d_refine_equivalence_classes.getSize();
         d_refine_equivalence_classes.resizeArray(eq_index+1);
      }

      if (data.d_oprefine.isNull()) {
         d_refine_equivalence_classes[eq_index].appendItem(data);
      } else {
         bool inserted = false;
         const int priority = data.d_oprefine->getOperatorPriority();
         typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator
            li(d_refine_equivalence_classes[eq_index]);
         while (!inserted && li) {
            if ((li().d_oprefine.isNull()) 
             || (li().d_oprefine->getOperatorPriority() >= priority)) {
               d_refine_equivalence_classes[eq_index].addItemBefore(li, data);
               inserted = true;
            }
            li++;
         }
         if (!inserted) {
            d_refine_equivalence_classes[eq_index].appendItem(data);
         }
      }

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

template<int DIM> bool RefineClasses<DIM>::checkRefineItem(
   const typename RefineClasses<DIM>::Data& data_item,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor) const
{

   bool item_good = true;

   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = descriptor;
   if (pd.isNull()) { 
      pd = hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor();
   }

   const int dst_id     = data_item.d_dst;
   const int src_id     = data_item.d_src;
   const int scratch_id = data_item.d_scratch;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "`Source' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (scratch_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "`Scratch' patch data id invalid (< 0!)" << std::endl);
   }

   tbox::Pointer< hier::PatchDataFactory<DIM> > dst_fact =
       pd->getPatchDataFactory(dst_id);
   tbox::Pointer< hier::PatchDataFactory<DIM> > src_fact =
       pd->getPatchDataFactory(src_id);
   tbox::Pointer< hier::PatchDataFactory<DIM> > scratch_fact =
       pd->getPatchDataFactory(scratch_id);

   if ( item_good && !(src_fact->validCopyTo(scratch_fact)) ) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "It is not a valid operation to copy from `Source' patch data \n"
                 << pd->mapIndexToName(src_id) << " to `Scratch' patch data "
                 << pd->mapIndexToName(scratch_id) << std::endl);
   }

   if ( item_good && !(scratch_fact->validCopyTo(dst_fact)) ) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "It is not a valid operation to copy from `Scratch' patch data \n"
                 << pd->mapIndexToName(scratch_id) << " to `Destination' patch data "
                 << pd->mapIndexToName(dst_fact) << std::endl);
   }

   const hier::IntVector<DIM>& scratch_gcw = scratch_fact->getGhostCellWidth();

   if (item_good && (dst_fact->getGhostCellWidth() > scratch_gcw)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                 << "`Destination' patch data " << pd->mapIndexToName(dst_id)
                 << " has a larger ghost cell width than \n"
                 << "`Scratch' patch data " << pd->mapIndexToName(scratch_id)
                 << "\n`Destination' ghost width = " 
                 << dst_fact->getGhostCellWidth()
                 << "\n`Scratch' ghost width = " << scratch_gcw << std::endl);
   }

   tbox::Pointer< RefineOperator<DIM> > refop = data_item.d_oprefine;
   if (item_good && !refop.isNull()) {
      if (refop->getStencilWidth() > scratch_gcw) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                    << "Refine operator "  << refop->getOperatorName()
                    << "\nhas larger stencil width than ghost cell width"
                    << "of `Scratch' patch data" << pd->mapIndexToName(scratch_id)
                    << "\noperator stencil width = " << refop->getStencilWidth()
                    << "\n`Scratch'  ghost width = " << scratch_gcw << std::endl);
      }
   }

   if (item_good && data_item.d_time_interpolate) {
      const int src_told_id = data_item.d_src_told; 
      const int src_tnew_id = data_item.d_src_tnew; 

      if (src_told_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                    << "`Source old' patch data id invalid (< 0!),\n"
                    << "yet a request has made to time interpolate" << std::endl);
      }
      if (item_good && src_tnew_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                    << "`Source new' patch data id invalid (< 0!),\n"
                    << "yet a request has made to time interpolate with them" 
                    << std::endl);
      }

      tbox::Pointer< hier::PatchDataFactory<DIM> > src_told_fact =
         pd->getPatchDataFactory(src_told_id);
      tbox::Pointer< hier::PatchDataFactory<DIM> > src_tnew_fact =
         pd->getPatchDataFactory(src_tnew_id);

      if ( item_good && typeid(*src_told_fact) != typeid(*src_fact) ) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                    << "`Source' patch data " << pd->mapIndexToName(src_id)
                    << " and `Source old' patch data " 
                    << pd->mapIndexToName(src_told_id)
                    << " have different patch data types, yet a request has"
                    << "\n been made to time interpolate with them" << std::endl);
      }

      if ( item_good && typeid(*src_tnew_fact) != typeid(*src_fact) ) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses<DIM>...\n"
                    << "`Source' patch data " << pd->mapIndexToName(src_id)
                    << " and `Source new' patch data "
                    << pd->mapIndexToName(src_tnew_id)
                    << " have different patch data types, yet a request has"
                    << "\n been made to time interpolate with them" << std::endl);
      }

   }

   return(item_good);

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

template<int DIM> bool RefineClasses<DIM>::checkConsistency(
   tbox::Pointer< RefineClasses<DIM> > test_classes,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor) const
{

   bool items_match = true; 

   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor();
   }

   if ( d_refine_equivalence_classes.getSize() !=
        test_classes->d_refine_equivalence_classes.getSize() ) {

      items_match = false;

   } else {

      int num_equiv_classes = d_refine_equivalence_classes.getSize();
      int eq_index = 0;
      while (items_match && eq_index < num_equiv_classes) {

         if ( d_refine_equivalence_classes[eq_index].getNumberItems() !=
              test_classes->
                 d_refine_equivalence_classes[eq_index].getNumberItems() ) {

            items_match = false;

         } else {

            typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator
               myli(d_refine_equivalence_classes[eq_index]); 
            typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator
               testli(test_classes->d_refine_equivalence_classes[eq_index]);
            while (items_match && myli) {

               items_match = checkPatchDataItemConsistency(
                             myli().d_dst, testli().d_dst, pd);

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                                myli().d_src, testli().d_src, pd);
               }

               if (items_match && myli().d_time_interpolate) {
                  items_match = checkPatchDataItemConsistency(
                                myli().d_src_told, testli().d_src_told, pd);
               }

               if (items_match && myli().d_time_interpolate) {
                  items_match = checkPatchDataItemConsistency(
                                myli().d_src_tnew, testli().d_src_tnew, pd);
               }

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                                myli().d_scratch, testli().d_scratch, pd);
               }

               if (items_match) {
                  items_match = (myli().d_fine_bdry_reps_var ==
                                 testli().d_fine_bdry_reps_var);
               }

               if (items_match) {
                  items_match = ( !myli().d_oprefine.isNull() ==
                                  !testli().d_oprefine.isNull() );
                  if (items_match && !myli().d_oprefine.isNull()) {
                     items_match = 
                        ( myli().d_oprefine->getStencilWidth() ==
                          testli().d_oprefine->getStencilWidth() );
                  }
               }

               if (items_match) {
                  items_match = ( myli().d_time_interpolate ==
                                  testli().d_time_interpolate);
                  if (items_match && myli().d_time_interpolate) {
                     items_match = ( typeid(*(myli().d_optime)) == 
                                     typeid(*(testli().d_optime)) );
                  }
               }

               myli++;
               testli++;

            } // while items in equivalence class match

         } // if number of items in equivalence class match

         eq_index++;

      } // while equivalence classes match

   } // else number of equivalence classes do not match 

   return(items_match);

}

/*
*************************************************************************
*                                                                       *
* Private member function to determine whether two patch data items     *
* are consistent.                                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool RefineClasses<DIM>::checkPatchDataItemConsistency(
   int item_id1,
   int item_id2,
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd) const
{

   bool items_match = ( (item_id1 >= 0) && (item_id2 >= 0) );

   if (items_match) {

      tbox::Pointer< hier::PatchDataFactory<DIM> > pdf1 =
         pd->getPatchDataFactory(item_id1);
      tbox::Pointer< hier::PatchDataFactory<DIM> > pdf2 =
         pd->getPatchDataFactory(item_id2);

      items_match = ( typeid(*pdf1) == typeid(*pdf2) );

      if (items_match) {
         items_match = ( pdf1->getGhostCellWidth() ==
                         pdf2->getGhostCellWidth() );
      }

   }

   return(items_match);

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

template<int DIM> int RefineClasses<DIM>::getEquivalenceClassIndex(
   const typename RefineClasses<DIM>::Data& data,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor) const
{

   int eq_index = -1;

   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor();
   }

   int dst_id = data.d_dst;
   int src_id = data.d_src;

   tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
      pd->getPatchDataFactory(dst_id);
   tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
      pd->getPatchDataFactory(src_id);

   hier::IntVector<DIM> dst_ghosts = dst_pdf->getGhostCellWidth();
   hier::IntVector<DIM> src_ghosts = src_pdf->getGhostCellWidth();

   int num_equiv_classes = d_refine_equivalence_classes.getSize();
   bool equiv_found = false;
   int nl = 0;
   while (!equiv_found && nl < num_equiv_classes) {

      bool dst_equiv = false;
      bool src_equiv = false;

      const typename RefineClasses<DIM>::Data& class_rep = getClassRepresentative(nl);

      int rep_dst_id = class_rep.d_dst;
      tbox::Pointer< hier::PatchDataFactory<DIM> > rep_dst_pdf =
         pd->getPatchDataFactory(rep_dst_id);
      hier::IntVector<DIM> rep_dst_ghosts =
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
         tbox::Pointer< hier::PatchDataFactory<DIM> > rep_src_pdf =
            pd->getPatchDataFactory(rep_src_id);
         hier::IntVector<DIM> rep_src_ghosts =
            rep_src_pdf->getGhostCellWidth();
         if ((src_ghosts == rep_src_ghosts) &&
             (typeid(*src_pdf) == typeid(*rep_src_pdf))) {
            src_equiv = true;
         }
      }

      /*
       * If destinations and sources are both equivalent, exit loop
       * and set return value to identify current equivalence class id.
       */
      if (dst_equiv && src_equiv) {
         eq_index = nl;
         equiv_found = true;
      }

      nl++;

   }

   return (eq_index);

}

/*
*************************************************************************
*									*
* Print the data in the refine item lists to the specified stream.      *
*									*
*************************************************************************
*/

template<int DIM> void RefineClasses<DIM>::printClassData(std::ostream& stream) const
{
   stream << "RefineClasses<DIM>::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < d_refine_equivalence_classes.getSize(); i++) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      for (typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator
           li(d_refine_equivalence_classes[i]); li; li++) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printRefineItem(stream, li());

         j++;
      }
      stream << std::endl;
   }

}

template<int DIM> void RefineClasses<DIM>::printRefineItem(
   std::ostream& stream,
   const typename RefineClasses<DIM>::Data& data) const
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
      stream << "time interpoate is false" << std::endl;
   } else {
      stream << "old source component:   " 
             << data.d_src_told << std::endl;
      stream << "new source component:   " 
             << data.d_src_tnew << std::endl;
      stream << "time interpolation operator name:          "
             << typeid(*data.d_optime).name()
             << std::endl;
   }
   stream << std::endl;
}

}
}
#endif
