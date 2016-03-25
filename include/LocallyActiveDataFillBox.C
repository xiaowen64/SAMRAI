//
// File:	LocallyActiveDataFillBox.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 481 $
// Modified:	$Date: 2005-07-21 13:50:43 -0700 (Thu, 21 Jul 2005) $
// Description:	Routines for "smart" boxlist ops in locally-active data comm ops
//

#ifndef included_xfer_LocallyActiveDataFillBox_C
#define included_xfer_LocallyActiveDataFillBox_C

#include "LocallyActiveDataFillBox.h"

#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
   namespace xfer {

template<int DIM>
LocallyActiveDataFillBox<DIM>::LocallyActiveDataFillBox(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data)
{
   d_box = &box;
   d_refine_var_data = var_data;
   d_refine_data = true;
}

template<int DIM>
LocallyActiveDataFillBox<DIM>::LocallyActiveDataFillBox(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data)
{
   d_box = &box;
   d_coarsen_var_data = var_data;
   d_refine_data = false;
}

template<int DIM>
LocallyActiveDataFillBox<DIM>::LocallyActiveDataFillBox(
   const LocallyActiveDataFillBox<DIM>& fill_box)
{
   d_box = fill_box.d_box;

   if (fill_box.d_refine_data) {
      d_refine_var_data = fill_box.d_refine_var_data;
      d_refine_data = true;
   } else {
      d_coarsen_var_data = fill_box.d_coarsen_var_data;
      d_refine_data = false;
   }
}

template<int DIM>
LocallyActiveDataFillBox<DIM>::~LocallyActiveDataFillBox()
{
   clearLocallyActiveFillBoxData();
}

template<int DIM>
const hier::Box<DIM>& LocallyActiveDataFillBox<DIM>::getBox() const
{
   return(*d_box);
}

template<int DIM>
const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& 
LocallyActiveDataFillBox<DIM>::getActiveRefineVarData() const
{
   if (!d_refine_data) {
      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::getActiveRefineVarData error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << endl);
   }
   return(d_refine_var_data);
}

template<int DIM>
const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>&
LocallyActiveDataFillBox<DIM>::getActiveCoarsenVarData() const
{
   if (d_refine_data) {
      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::getActiveCoarsenVarData error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << endl);
   }
   return(d_coarsen_var_data);
}

template<int DIM>
void LocallyActiveDataFillBox<DIM>::clearLocallyActiveFillBoxData()
{
   d_box = (hier::Box<DIM>*)NULL;

   for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
      lavdi(d_refine_var_data); lavdi; lavdi++) {
      lavdi() = (const typename xfer::RefineClasses<DIM>::Data*)NULL;
   }
   d_refine_var_data.clearItems();

   for (typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator
      lavdi(d_coarsen_var_data); lavdi; lavdi++) {
      lavdi() = (const typename xfer::CoarsenClasses<DIM>::Data*)NULL;
   }
   d_coarsen_var_data.clearItems();
}

template<int DIM>
void LocallyActiveDataFillBox<DIM>::printClassData(ostream& os) const
{
   const hier::Box<DIM> print_box(*d_box);
   os << "\n Box = " << print_box << endl;
   if (d_refine_data) {
      if (d_refine_var_data.getNumberItems() == 0) {
         os << "  no refine var ids ";
      } else {
         os << "  refine var ids = ...";

         for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
              lavdi(d_refine_var_data); lavdi; lavdi++) {
            os << "\n     dst, src, scratch = " 
               << lavdi()->d_dst << " , " << lavdi()->d_src << " , " << lavdi()->d_scratch;
         }
      }
   } else {
      if (d_coarsen_var_data.getNumberItems() == 0) {
         os << "  no coarsen var ids ";
      } else {
         os << "  coarsen var ids = ...";

         for (typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator
              lavdi(d_coarsen_var_data); lavdi; lavdi++) {
            os << "\n     dst, src = "
               << lavdi()->d_dst << " , " << lavdi()->d_src;
         }
      }
   }
}

template<int DIM>
bool LocallyActiveDataFillBox<DIM>::checkData(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data,
   ostream& os) const
{
   bool ret_val = true;

   if (!d_refine_data) {

      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::checkData error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << endl);
      ret_val = false;

   } else {

      if (box != *d_box) {
         ret_val = false;
         os << "LocallyActiveDataFillBox<DIM>::checkData..."
            << "   this box = " << *d_box
            << "\n   input box = " << box << endl;
      }

      if (ret_val) {
         if (var_data.getNumberItems() != d_refine_var_data.getNumberItems()) {
            ret_val = false;
            os << "LocallyActiveDataFillBox<DIM>::checkData..."
               << "   # var items in this set = " << d_refine_var_data.getNumberItems()
               << "\n   # var items in input = " << var_data.getNumberItems() << endl;
         } else {
            typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator 
               lavdi_input(var_data);
            typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator 
               lavdi(d_refine_var_data);
            int cnt = 0;
            for ( ; lavdi; lavdi++) {
               if ( (lavdi_input()->d_dst != lavdi()->d_dst) ||
                    (lavdi_input()->d_src != lavdi()->d_src) ||
                    (lavdi_input()->d_scratch != lavdi()->d_scratch) ) {
                  ret_val = false;
                  os << "LocallyActiveDataFillBox<DIM>::checkData..."
                     << "   list item # " << cnt
                     << "\n    this var ids(dst,src,scr) = "
                        << lavdi()->d_dst << " , " 
                        << lavdi()->d_src << " , " 
                        << lavdi()->d_scratch
                     << "\n    input var ids(dst,src,scr) = "
                        << lavdi_input()->d_dst << " , " 
                        << lavdi_input()->d_src << " , " 
                        << lavdi_input()->d_scratch;
               }
               lavdi_input++;
               cnt++;
            }
         }
      }

   }

   return(ret_val);
}

template<int DIM>
bool LocallyActiveDataFillBox<DIM>::checkData(
   const hier::Box<DIM>& box,
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data,
   ostream& os) const
{
   bool ret_val = true;

   if (d_refine_data) {

      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::checkData error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << endl);
      ret_val = false;

   } else {

      if (box != *d_box) {
         ret_val = false;
         os << "LocallyActiveDataFillBox<DIM>::checkData..."
            << "   this box = " << *d_box
            << "\n   input box = " << box << endl;
      }

      if (ret_val) {
         if (var_data.getNumberItems() != d_coarsen_var_data.getNumberItems()) {
            ret_val = false;
            os << "LocallyActiveDataFillBox<DIM>::checkData..."
               << "   # var items in this set = " << d_coarsen_var_data.getNumberItems()
               << "\n   # var items in input = " << var_data.getNumberItems() << endl;
         } else {
            typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
               lavdi_input(var_data);
            typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
               lavdi(d_coarsen_var_data);
            int cnt = 0;
            for ( ; lavdi; lavdi++) {
               if ( (lavdi_input()->d_dst != lavdi()->d_dst) ||
                    (lavdi_input()->d_src != lavdi()->d_src) ) {
                  ret_val = false;
                  os << "LocallyActiveDataFillBox<DIM>::checkData..."
                     << "   list item # " << cnt
                     << "\n    this var ids(dst,src) = "
                        << lavdi()->d_dst << " , " 
                        << lavdi()->d_src 
                     << "\n    input var ids(dst,src) = "
                        << lavdi_input()->d_dst << " , " 
                        << lavdi_input()->d_src;
               }
               lavdi_input++;
               cnt++;
            }
         }
      }

   }

   return(ret_val);
}

}
}
#endif
