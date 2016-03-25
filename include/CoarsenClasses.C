//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenClasses.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2043 $
// Modified:	$LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description:	Simple structure for managing coarsening data in equivalence classes.
//

#ifndef included_xfer_CoarsenClasses_C
#define included_xfer_CoarsenClasses_C

#include <typeinfo>

#include "CoarsenClasses.h"

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
* Constructor sets boolean for filling coarse data and creates new      *
* array of equivalence classes.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CoarsenClasses<DIM>::CoarsenClasses(bool fill_coarse_data)
{
   d_fill_coarse_data = fill_coarse_data;
   d_coarsen_equivalence_classes.resizeArray(0);
}

/*
*************************************************************************
*									*
* The destructor implicitly deletes the item storage associated with	*
* the equivalence classes (and also the coarsen algorithm).		*
*									*
*************************************************************************
*/

template<int DIM>  CoarsenClasses<DIM>::~CoarsenClasses()
{
   d_coarsen_equivalence_classes.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Return number of equivalence classes.                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> int CoarsenClasses<DIM>::getNumberOfEquivalenceClasses() const
{
   return(d_coarsen_equivalence_classes.getSize());
}

/*
*************************************************************************
*                                                                       *
* Return representative item for given equivalence class (first in list)*
*                                                                       *
*************************************************************************
*/

template<int DIM> const typename CoarsenClasses<DIM>::Data& 
   CoarsenClasses<DIM>::getClassRepresentative(int equiv_class_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (equiv_class_id >= 0) && 
           (equiv_class_id < d_coarsen_equivalence_classes.getSize()) );
#endif 
   return(d_coarsen_equivalence_classes[equiv_class_id].getFirstItem());
}

/*
*************************************************************************
*                                                                       *
* Return iterator for list of coaren items for given equivalence class  *
*                                                                       *
*************************************************************************
*/

template<int DIM> typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator 
CoarsenClasses<DIM>::getIterator(int equiv_class_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (equiv_class_id >= 0) && 
           (equiv_class_id < d_coarsen_equivalence_classes.getSize()) );
#endif 
   return(typename tbox::List<typename CoarsenClasses<DIM>::Data>::
          Iterator(d_coarsen_equivalence_classes[equiv_class_id]));
}

/*
*************************************************************************
*									*
* Insert a data item into the coarsen data list for the proper          *
* equivalence class in sorted order by ascending operator priority.     *
*									*
*************************************************************************
*/

template<int DIM> void CoarsenClasses<DIM>::insertEquivalenceClassItem(
   const typename CoarsenClasses<DIM>::Data& data,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor)
{

   if (!checkCoarsenItem(data, descriptor)) {
      tbox::perr << "Bad refine class data passed to "
           << "CoarsenClasses<DIM>::insertEquivalenceClassItem\n";
      printCoarsenItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = d_coarsen_equivalence_classes.getSize();
         d_coarsen_equivalence_classes.resizeArray(eq_index+1);
      }

      if (data.d_opcoarsen.isNull()) {
         d_coarsen_equivalence_classes[eq_index].appendItem(data);
      } else {
         bool inserted = false;
         const int priority = data.d_opcoarsen->getOperatorPriority();
         typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator
            li(d_coarsen_equivalence_classes[eq_index]);
         while (!inserted && li) {
            if ( (li().d_opcoarsen.isNull()) || 
                 (li().d_opcoarsen->getOperatorPriority() >= priority) ) {
               d_coarsen_equivalence_classes[eq_index].addItemBefore(li, data);
               inserted = true;
            }
            li++;
         }
         if (!inserted) {
            d_coarsen_equivalence_classes[eq_index].appendItem(data);
         }
      }

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

template<int DIM> bool CoarsenClasses<DIM>::checkCoarsenItem(
   const typename CoarsenClasses<DIM>::Data& data_item,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor) const
{

   bool item_good = true;

   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = descriptor;
   if (pd.isNull()) { 
      pd = hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor();
   }

   const int dst_id = data_item.d_dst;
   const int src_id = data_item.d_src;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses<DIM>...\n"
                 << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses<DIM>...\n"
                 << "`Source' patch data id invalid (< 0!)" << std::endl);
   }

   tbox::Pointer< hier::PatchDataFactory<DIM> > dfact =
      pd->getPatchDataFactory(dst_id);
   tbox::Pointer< hier::PatchDataFactory<DIM> > sfact =
      pd->getPatchDataFactory(src_id);

   if ( item_good && !(sfact->validCopyTo(dfact)) ) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses<DIM>...\n"
                 << "It is not a valid operation to copy from `Source' patch data \n" 
                 << pd->mapIndexToName(src_id) << " to `Destination' patch data " 
                 << pd->mapIndexToName(dst_id) << std::endl);
   }

   tbox::Pointer< CoarsenOperator<DIM> > coarsop = data_item.d_opcoarsen;
   if (item_good && !coarsop.isNull()) {
      if (coarsop->getStencilWidth() > sfact->getGhostCellWidth()) {
         item_good = false;
         TBOX_ERROR("Bad data given to CoarsenClasses<DIM>...\n"
                    << "Coarsen operator "  << coarsop->getOperatorName()
                    << "\nhas larger stencil width than ghost cell width"
                    << "of `Source' patch data" << pd->mapIndexToName(src_id)
                    << "\noperator stencil width = " << coarsop->getStencilWidth()
                    << "\n`Source'  ghost width = " 
                    << sfact->getGhostCellWidth()
                    << std::endl);
      }
   }

   return(item_good);

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

template<int DIM> bool CoarsenClasses<DIM>::checkConsistency(
   tbox::Pointer< CoarsenClasses<DIM> > test_classes,
   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor) const
{

   bool items_match = true;

   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = descriptor;
   if (pd.isNull()) {
      pd = hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor();
   }

   if ( d_coarsen_equivalence_classes.getSize() !=
        test_classes->d_coarsen_equivalence_classes.getSize() ) {

      items_match = false;

   } else {

      int num_equiv_classes = d_coarsen_equivalence_classes.getSize();
      int eq_index = 0;
      while (items_match && eq_index < num_equiv_classes) {

         if ( d_coarsen_equivalence_classes[eq_index].getNumberOfItems() !=
              test_classes->
                 d_coarsen_equivalence_classes[eq_index].getNumberOfItems() ) {

            items_match = false;

         } else {

            typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator
               myli(d_coarsen_equivalence_classes[eq_index]);
            typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator
               testli(test_classes->d_coarsen_equivalence_classes[eq_index]);
            while (items_match && myli) {

               items_match = checkPatchDataItemConsistency(
                             myli().d_dst, testli().d_dst, pd);

               if (items_match) {
                  items_match = checkPatchDataItemConsistency(
                                myli().d_src, testli().d_src, pd);
               }

               if (items_match) {
                  items_match = (myli().d_fine_bdry_reps_var ==
                                 testli().d_fine_bdry_reps_var);
               }

               if (items_match) {
                  items_match = (myli().d_gcw_to_coarsen ==
                                 testli().d_gcw_to_coarsen);
               }

               if (items_match) {
                  items_match = ( !myli().d_opcoarsen.isNull() ==
                                  !testli().d_opcoarsen.isNull() );
                  if (items_match && !myli().d_opcoarsen.isNull()) {
                     items_match =
                        ( myli().d_opcoarsen->getStencilWidth() ==
                          testli().d_opcoarsen->getStencilWidth() );
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

template<int DIM> bool CoarsenClasses<DIM>::checkPatchDataItemConsistency(
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

template<int DIM> int CoarsenClasses<DIM>::getEquivalenceClassIndex(
   const typename CoarsenClasses<DIM>::Data& data,
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

   int num_equiv_classes = d_coarsen_equivalence_classes.getSize();
   bool equiv_found = false;
   int nl = 0;
   while (!equiv_found && nl < num_equiv_classes) {

      bool dst_equiv = false;
      bool src_equiv = false;

      const typename CoarsenClasses<DIM>::Data& class_rep = getClassRepresentative(nl);

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
*                                                                       *
* Return the number of items in the specified equivalence class.        *
*                                                                       *
*************************************************************************
*/

template<int DIM> int CoarsenClasses<DIM>::getNumberOfItemsInEquivalenceClass(
   int equiv_class_id) const
{
   return (d_coarsen_equivalence_classes[equiv_class_id].size());
}

/*
*************************************************************************
*									*
* Print the data in the coarsen item lists to the specified stream.     *
*									*
*************************************************************************
*/

template<int DIM> void CoarsenClasses<DIM>::printClassData(std::ostream& stream) const
{
   stream << "CoarsenClasses<DIM>::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < d_coarsen_equivalence_classes.getSize(); i++) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      for (typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator
           li(d_coarsen_equivalence_classes[i]); li; li++) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printCoarsenItem(stream, li());

         j++;
      }
      stream << std::endl;
   }

}

template<int DIM> void CoarsenClasses<DIM>::printCoarsenItem(
   std::ostream& stream,
   const typename CoarsenClasses<DIM>::Data& data) const
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
