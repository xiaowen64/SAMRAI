//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/templates/special/tbox_List-xfer_CoarsenClasses_NDIM__Data-PtrX.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	special template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "CoarsenClasses.h"

namespace SAMRAI {
        namespace tbox {

#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< const xfer::CoarsenClasses<NDIM>::Data* > *tbox::ListNode< const xfer::CoarsenClasses<NDIM>::Data* >::s_free_list=0;
bool tbox::ListNode< const xfer::CoarsenClasses<NDIM>::Data* >::s_registered_callback=false;
#endif

template class tbox::List< const xfer::CoarsenClasses<NDIM>::Data* >;
template class tbox::ListIterator< const xfer::CoarsenClasses<NDIM>::Data* >;
template class tbox::ListNode< const xfer::CoarsenClasses<NDIM>::Data* >;

}
}
