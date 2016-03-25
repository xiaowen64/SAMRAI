//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataFillBoxSet.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for "smart" boxlist ops in locally-active comm schedules
//

#ifndef included_xfer_LocallyActiveDataFillBoxSet_C
#define included_xfer_LocallyActiveDataFillBoxSet_C

#include "LocallyActiveDataFillBoxSet.h"

#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

namespace SAMRAI {
   namespace xfer {

//#define LADFBS_CHECK_UNION
#undef LADFBS_CHECK_UNION

template<int DIM>
LocallyActiveDataFillBoxSet<DIM>::LocallyActiveDataFillBoxSet(
   const LocallyActiveDataFillBoxSet<DIM>& fill_box_set)
: xfer::FillBoxSet<DIM>(fill_box_set)
{
   typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
      ladfbi(fill_box_set.d_locally_active_boxes);

   const hier::BoxList<DIM>& boxes = FillBoxSet<DIM>::getBoxList();

   if (fill_box_set.d_refine_data) {
      d_refine_data = true;
      for (typename hier::BoxList<DIM>::Iterator lb(boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveRefineVarData());
         d_locally_active_boxes.appendItem(fill_box);

         ladfbi++;
      }
      d_union_refine_var_data = fill_box_set.getUnionActiveRefineVarData();
   } else {
      d_refine_data = false;
      for (typename hier::BoxList<DIM>::Iterator lb(boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveCoarsenVarData());
         d_locally_active_boxes.appendItem(fill_box);

         ladfbi++;
      }
      d_union_coarsen_var_data = fill_box_set.getUnionActiveCoarsenVarData();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION 
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
LocallyActiveDataFillBoxSet<DIM>::LocallyActiveDataFillBoxSet()
: xfer::FillBoxSet<DIM>()
{
   clearLocallyActiveFillBoxData();
   d_refine_data = false;
}

template<int DIM>
LocallyActiveDataFillBoxSet<DIM>::~LocallyActiveDataFillBoxSet()
{
   clearLocallyActiveFillBoxData();
}

template<int DIM>
int LocallyActiveDataFillBoxSet<DIM>::getNumberOfBoxes() const
{
   return(d_locally_active_boxes.getNumberItems());
}

template<int DIM>
const tbox::List< xfer::LocallyActiveDataFillBox<DIM> >& 
LocallyActiveDataFillBoxSet<DIM>::getLocallyActiveDataBoxes() const
{
   return(d_locally_active_boxes);
}

template<int DIM>
const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>&
LocallyActiveDataFillBoxSet<DIM>::getUnionActiveRefineVarData() const
{
   if (!d_refine_data) {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::getUnionActiveRefineVarData error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
   return(d_union_refine_var_data);
}

template<int DIM>
const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>&
LocallyActiveDataFillBoxSet<DIM>::getUnionActiveCoarsenVarData() const
{
   if (d_refine_data) {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::getActiveCoarsenVarData error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << std::endl);
   }
   return(d_union_coarsen_var_data);
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::setTo(
   const LocallyActiveDataFillBoxSet<DIM>& la_fill_box_set)
{
   if (this != &la_fill_box_set) {

      this->clearLocallyActiveFillBoxData();

      resetFillBoxes(la_fill_box_set.getBoxList());
      const hier::BoxList<DIM>& new_boxes = FillBoxSet<DIM>::getBoxList();

      typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
         ladfbi(la_fill_box_set.d_locally_active_boxes);

      if (la_fill_box_set.d_refine_data) {
         d_refine_data = true;
         for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
            xfer::LocallyActiveDataFillBox<DIM> 
               fill_box(lb(), ladfbi().getActiveRefineVarData());
            d_locally_active_boxes.appendItem(fill_box);
            ladfbi++;
         }
         d_union_refine_var_data = la_fill_box_set.getUnionActiveRefineVarData();
      } else {
         d_refine_data = false;
         for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
            xfer::LocallyActiveDataFillBox<DIM> 
               fill_box(lb(), ladfbi().getActiveCoarsenVarData());
            d_locally_active_boxes.appendItem(fill_box);
            ladfbi++;
         }
         d_union_coarsen_var_data = la_fill_box_set.getUnionActiveCoarsenVarData();
      }

   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::resetLocallyActiveFillBoxes(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data)
{
   if (d_refine_data || d_locally_active_boxes.getNumberItems() == 0) {

      d_refine_data = true;

      resetFillBoxes(box);
      clearLocallyActiveFillBoxData(); 

      xfer::LocallyActiveDataFillBox<DIM> fill_box(FillBoxSet<DIM>::getBoundingBox(),
                                                   var_data);

      d_locally_active_boxes.appendItem(fill_box);

      d_union_refine_var_data = var_data;

   } else {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::resetLocallyActiveFillBoxes error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::resetLocallyActiveFillBoxes(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data)
{
   if (!d_refine_data || d_locally_active_boxes.getNumberItems() == 0) {

      d_refine_data = false;

      resetFillBoxes(box);
      clearLocallyActiveFillBoxData();

      xfer::LocallyActiveDataFillBox<DIM> fill_box(FillBoxSet<DIM>::getBoundingBox(),
                                                   var_data);

      d_locally_active_boxes.appendItem(fill_box);

      d_union_coarsen_var_data = var_data;

   } else {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::resetLocallyActiveFillBoxes error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::addLocallyActiveFillBox(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data)
{
   if (d_refine_data || d_locally_active_boxes.getNumberItems() == 0) {
 
      d_refine_data = true;

      xfer::FillBoxSet<DIM>::addFillBox(box);

      xfer::LocallyActiveDataFillBox<DIM> fill_box(FillBoxSet<DIM>::getBoxList().getFirstItem(),
                                                   var_data);

      d_locally_active_boxes.addItem(fill_box);

      tbox::List<const typename xfer::RefineClasses<DIM>::Data*> tmp_list;
      mergeLists(tmp_list, d_union_refine_var_data, var_data);
      d_union_refine_var_data = tmp_list;

   } else {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::addFillBox error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::addLocallyActiveFillBox(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data)
{
   if (!d_refine_data || d_locally_active_boxes.getNumberItems() == 0) {

      d_refine_data = false;

      xfer::FillBoxSet<DIM>::addFillBox(box);

      xfer::LocallyActiveDataFillBox<DIM> fill_box(FillBoxSet<DIM>::getBoxList().getFirstItem(),
                                                   var_data);

      d_locally_active_boxes.addItem(fill_box);

      tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> tmp_list;
      mergeLists(tmp_list, d_union_coarsen_var_data, var_data);
      d_union_coarsen_var_data = tmp_list;

   } else {
      TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::addLocallyActiveFillBox error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::intersectBoxes(const hier::Box<DIM>& box)
{
   hier::BoxList<DIM> intersection_boxes;
   tbox::List< xfer::LocallyActiveDataFillBox<DIM> > locally_active_intersection_boxes;

   while ( !(d_locally_active_boxes.isEmpty()) ) {
      xfer::LocallyActiveDataFillBox<DIM> tryme = d_locally_active_boxes.getFirstItem();
      d_locally_active_boxes.removeFirstItem();
      hier::Box<DIM> intersection = tryme.getBox() * box;
      if (!intersection.empty()) {
         intersection_boxes.appendItem(intersection);
         if (d_refine_data) {
            xfer::LocallyActiveDataFillBox<DIM> fill_box(intersection,
                                                         tryme.getActiveRefineVarData());
            locally_active_intersection_boxes.appendItem(fill_box);
         } else {
            xfer::LocallyActiveDataFillBox<DIM> fill_box(intersection,
                                                         tryme.getActiveCoarsenVarData());
            locally_active_intersection_boxes.appendItem(fill_box);
         }
      }
   }

   clearLocallyActiveFillBoxData();
   resetFillBoxes(intersection_boxes);

   const hier::BoxList<DIM>& new_boxes = FillBoxSet<DIM>::getBoxList();
   typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
      ladfbi(locally_active_intersection_boxes);

   if (d_refine_data) {
      for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveRefineVarData());
         d_locally_active_boxes.appendItem(fill_box);

         tbox::List<const typename xfer::RefineClasses<DIM>::Data*> tmp_list;
         mergeLists(tmp_list, d_union_refine_var_data, ladfbi().getActiveRefineVarData());
         d_union_refine_var_data = tmp_list;

         ladfbi++;
      }
   } else {
      for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveCoarsenVarData());
         d_locally_active_boxes.appendItem(fill_box);

         tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> tmp_list;
         mergeLists(tmp_list, d_union_coarsen_var_data, ladfbi().getActiveCoarsenVarData());
         d_union_coarsen_var_data = tmp_list;

         ladfbi++;
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::intersectBoxes(const hier::BoxList<DIM>& boxes)
{
   hier::BoxList<DIM> intersection_boxes;
   tbox::List< xfer::LocallyActiveDataFillBox<DIM> > locally_active_intersection_boxes; 

   while ( !(d_locally_active_boxes.isEmpty()) ) {
      xfer::LocallyActiveDataFillBox<DIM> tryme = d_locally_active_boxes.getFirstItem();
      d_locally_active_boxes.removeFirstItem();
      for (typename hier::BoxList<DIM>::Iterator i(boxes); i; i++) {
         hier::Box<DIM> intersection = tryme.getBox() * i();
         if (!intersection.empty()) {
            intersection_boxes.appendItem(intersection);
            if (d_refine_data) {
               xfer::LocallyActiveDataFillBox<DIM> 
                  fill_box(intersection, tryme.getActiveRefineVarData());
               locally_active_intersection_boxes.appendItem(fill_box);
            } else {
               xfer::LocallyActiveDataFillBox<DIM> 
                  fill_box(intersection, tryme.getActiveCoarsenVarData());
               locally_active_intersection_boxes.appendItem(fill_box);
            }
         }
      }
   }

   clearLocallyActiveFillBoxData();
   resetFillBoxes(intersection_boxes);

   const hier::BoxList<DIM>& new_boxes = FillBoxSet<DIM>::getBoxList();
   typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
      ladfbi(locally_active_intersection_boxes);

   if (d_refine_data) {
      for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveRefineVarData());
         d_locally_active_boxes.appendItem(fill_box);

         tbox::List<const typename xfer::RefineClasses<DIM>::Data*> tmp_list;
         mergeLists(tmp_list, d_union_refine_var_data, ladfbi().getActiveRefineVarData());
         d_union_refine_var_data = tmp_list;

         ladfbi++;
      }   
   } else {
      for (typename hier::BoxList<DIM>::Iterator lb(new_boxes); lb; lb++) {
         xfer::LocallyActiveDataFillBox<DIM> 
            fill_box(lb(), ladfbi().getActiveCoarsenVarData());
         d_locally_active_boxes.appendItem(fill_box);

         tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> tmp_list;
         mergeLists(tmp_list, d_union_coarsen_var_data, ladfbi().getActiveCoarsenVarData());
         d_union_coarsen_var_data = tmp_list;

         ladfbi++;
      }   
   }
#ifdef DEBUG_CHECK_ASSERTIONS
#ifdef LADFBS_CHECK_UNION
   (void) checkUnion(tbox::plog);
#endif
#endif
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::printClassData(std::ostream& os) const
{
   for (typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
        ladfbi(d_locally_active_boxes); ladfbi; ladfbi++) {
      os << std::endl;
      ladfbi().printClassData(os);
   }
   os << std::endl;
   os << "\nxfer::FillBoxSet<DIM> object data...." << std::endl;
   xfer::FillBoxSet<DIM>::print(os);

   if (d_refine_data) {
      os << "\nd_union_refine_var_data = ..." << std::endl;
      printRefineVarListItems(d_union_refine_var_data, os);
   } else {
      os << "\nd_union_coarsen_var_data = ..." << std::endl;
      if (d_union_coarsen_var_data.getNumberItems() == 0) {
         os << "  no coarsen var ids in set ";
      } else {
         os << "  coarsen var ids = ...";

         for (typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator
              lavdi(d_union_coarsen_var_data); lavdi; lavdi++) {
            os << "\n     dst, src = "
               << lavdi()->d_dst << " , " << lavdi()->d_src;
         }
      }
   }

   os << std::endl;
   
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::clearLocallyActiveFillBoxData()
{
   for (typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator
        ladfbi(d_locally_active_boxes); ladfbi; ladfbi++) {
      ladfbi().clearLocallyActiveFillBoxData();
   }
   d_locally_active_boxes.clearItems();

   d_union_refine_var_data.clearItems();
   d_union_coarsen_var_data.clearItems();

}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::mergeLists(
   tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& outlist, 
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& inlist_a,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& inlist_b) const
{
   outlist.clearItems();

   typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator 
      ila(inlist_a);
   typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator 
      ilb(inlist_b);

   int last_add = tbox::MathUtilities<int>::getMin();

   while (ila && ilb) {
      int value_a = ila()->d_tag;
      int value_b = ilb()->d_tag;
      int value =
         tbox::MathUtilities<int>::Max( last_add,
                                        tbox::MathUtilities<int>::Min( value_a, 
                                                                       value_b) );
       
      if (value > last_add) {
         if (value_a < value_b) {
            outlist.appendItem(ila());
         } else {
            outlist.appendItem(ilb());
         }
         last_add = value;
      }

      while ( ila && (ila()->d_tag <= last_add) ) ila++;
      while ( ilb && (ilb()->d_tag <= last_add) ) ilb++;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !(ila && ilb) );
#endif

   if (!outlist.isEmpty()) {
      last_add = outlist.getLastItem()->d_tag;
   }

   typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator  
      ilremain( (ila ? ila : ilb) );
   while (ilremain) {
      int value = ilremain()->d_tag;
      if (value > last_add) {
         outlist.appendItem(ilremain());
         last_add = value;
      }
      ilremain++;
   }

}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::mergeLists(
   tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& outlist, 
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& inlist_a,
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& inlist_b) const
{
   outlist.clearItems();

   typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
      ila(inlist_a);
   typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
      ilb(inlist_b);

   int last_add = tbox::MathUtilities<int>::getMin();

   while (ila && ilb) {
      int value_a = ila()->d_tag;
      int value_b = ilb()->d_tag;
      int value =
         tbox::MathUtilities<int>::Max( last_add,
                                        tbox::MathUtilities<int>::Min(value_a, 
                                                                      value_b) );

      if (value > last_add) {
         if (value_a < value_b) {
            outlist.appendItem(ila());
         } else {
            outlist.appendItem(ilb());
         }
         last_add = value;
      }

      while ( ila && (ila()->d_tag <= last_add) ) ila++;
      while ( ilb && (ilb()->d_tag <= last_add) ) ilb++;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !(ila && ilb) );
#endif

   if (!outlist.isEmpty()) {
      last_add = outlist.getLastItem()->d_tag;
   }

   typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
      ilremain( (ila ? ila : ilb) );
   while (ilremain) {
      int value = ilremain()->d_tag;
      if (value > last_add) {
         outlist.appendItem(ilremain());
         last_add = value;
      }
      ilremain++;
   }

}

template<int DIM>
bool LocallyActiveDataFillBoxSet<DIM>::checkUnion(std::ostream& os) const
{ 
   bool union_match = true;

   tbox::List<const typename xfer::RefineClasses<DIM>::Data*> test_union;

   if (!d_refine_data) {

       TBOX_ERROR("LocallyActiveDataFillBoxSet<DIM>::checkUnion() error... " 
                  << "not set up for refine data." << std::endl);

   } else {

      typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator
         ladfbi(d_locally_active_boxes);

      os << "\n In checkUnion -- Computing union..." << std::endl;
      int fb_count = 0;
      for ( ; ladfbi; ladfbi++) {
         os << "\n ladfbi = " << fb_count << " : fill box active data = ..." << std::endl;
         printRefineVarListItems(ladfbi().getActiveRefineVarData(), os);
         os << std::endl;

         tbox::List<const typename xfer::RefineClasses<DIM>::Data*> tmp_list;
         mergeLists(tmp_list, test_union, ladfbi().getActiveRefineVarData());
         test_union = tmp_list;

         os << "\n test_union active data = ..." << std::endl;
         printRefineVarListItems(test_union, os);
         os << std::endl;
         fb_count++;
      }

      if (test_union.size() != d_union_refine_var_data.size()) {
         union_match = false;
      } else {

         typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
            il_test(test_union);
         typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
            il_union(d_union_refine_var_data); 
         while (il_test && il_union && union_match) {
            if (il_test()->d_tag != il_union()->d_tag) {
               union_match = false;
            }
            il_test++;
            il_union++;
         }

      }

   }

   if (!union_match) {
      os << "\n\n UNION DOES NOT MATCH!" << std::endl;
      os << "   Here is fill box set..." << std::endl;
      printClassData(os);
      os << "\n   Here is computed union..." << std::endl;
      printRefineVarListItems(test_union, os);
      os << std::endl;
      TBOX_ERROR("UNION DOES NOT MATCH!" << std::endl);
   }

   return(union_match);
}

template<int DIM>
void LocallyActiveDataFillBoxSet<DIM>::printRefineVarListItems(
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& inlist,
   std::ostream& os) const
{
   os << "list size = " << inlist.size() << std::endl;
   if (inlist.size() == 0) {
      os << "  no refine var ids in set ";
   } else {
      os << "  refine var ids = ...";

      for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
           lavdi(inlist); lavdi; lavdi++) {
         os << "\n     d_tag: dst, src, scratch = "
            << lavdi()->d_tag << " : " 
            << lavdi()->d_dst << " , " << lavdi()->d_src << " , " << lavdi()->d_scratch;
      }
   }
   os << std::endl;
}

}
}

#endif
