//
// File:	tbox_List-xfer_CoarsenClasses_NDIM__Data-PtrX.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 302 $
// Modified:	$Date: 2005-04-25 10:35:34 -0700 (Mon, 25 Apr 2005) $
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
