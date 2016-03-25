//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataFillBox.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for "smart" boxlist ops in locally-active data comm ops
//

#ifndef included_xfer_LocallyActiveDataFillBox_C
#define included_xfer_LocallyActiveDataFillBox_C

#include "LocallyActiveDataFillBox.h"

#include "tbox/Utilities.h"

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
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << std::endl);
   }
   return(d_refine_var_data);
}

template<int DIM>
const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>&
LocallyActiveDataFillBox<DIM>::getActiveCoarsenVarData() const
{
   if (d_refine_data) {
      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::getActiveCoarsenVarData error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << std::endl);
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
void LocallyActiveDataFillBox<DIM>::printClassData(std::ostream& os) const
{
   const hier::Box<DIM> print_box(*d_box);
   os << "\n Box = " << print_box << std::endl;
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
   std::ostream& os) const
{
   bool ret_val = true;

   if (!d_refine_data) {

      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::checkData error... "
                  << "\nobject is setup for xfer::CoarsenClasses<DIM>::Data" << std::endl);
      ret_val = false;

   } else {

      if (box != *d_box) {
         ret_val = false;
         os << "LocallyActiveDataFillBox<DIM>::checkData..."
            << "   this box = " << *d_box
            << "\n   input box = " << box << std::endl;
      }

      if (ret_val) {
         if (var_data.getNumberItems() != d_refine_var_data.getNumberItems()) {
            ret_val = false;
            os << "LocallyActiveDataFillBox<DIM>::checkData..."
               << "   # var items in this set = " << d_refine_var_data.getNumberItems()
               << "\n   # var items in input = " << var_data.getNumberItems() << std::endl;
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
   std::ostream& os) const
{
   bool ret_val = true;

   if (d_refine_data) {

      TBOX_ERROR("LocallyActiveDataFillBox<DIM>::checkData error... "
                  << "\nobject is setup for xfer::RefineClasses<DIM>::Data" << std::endl);
      ret_val = false;

   } else {

      if (box != *d_box) {
         ret_val = false;
         os << "LocallyActiveDataFillBox<DIM>::checkData..."
            << "   this box = " << *d_box
            << "\n   input box = " << box << std::endl;
      }

      if (ret_val) {
         if (var_data.getNumberItems() != d_coarsen_var_data.getNumberItems()) {
            ret_val = false;
            os << "LocallyActiveDataFillBox<DIM>::checkData..."
               << "   # var items in this set = " << d_coarsen_var_data.getNumberItems()
               << "\n   # var items in input = " << var_data.getNumberItems() << std::endl;
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
