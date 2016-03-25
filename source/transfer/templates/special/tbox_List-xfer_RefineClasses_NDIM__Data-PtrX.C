//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/templates/special/tbox_List-xfer_RefineClasses_NDIM__Data-PtrX.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	special template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "RefineClasses.h"

namespace SAMRAI {
        namespace tbox {

#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< const xfer::RefineClasses<NDIM>::Data* > *tbox::ListNode< const xfer::RefineClasses<NDIM>::Data* >::s_free_list=0;
bool tbox::ListNode< const xfer::RefineClasses<NDIM>::Data* >::s_registered_callback=false;
#endif

template class tbox::List< const xfer::RefineClasses<NDIM>::Data* >;
template class tbox::ListIterator< const xfer::RefineClasses<NDIM>::Data* >;
template class tbox::ListNode< const xfer::RefineClasses<NDIM>::Data* >;

}
}
