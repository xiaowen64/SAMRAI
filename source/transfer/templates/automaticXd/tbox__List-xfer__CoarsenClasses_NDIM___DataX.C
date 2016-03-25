//
// File:	tbox__List-xfer__CoarsenClasses_NDIM___DataX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "CoarsenClasses.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< xfer::CoarsenClasses<NDIM>::Data > *tbox::ListNode< xfer::CoarsenClasses<NDIM>::Data >::s_free_list=0;
bool tbox::ListNode< xfer::CoarsenClasses<NDIM>::Data >::s_registered_callback=false;
#endif
template class tbox::List< xfer::CoarsenClasses<NDIM>::Data >;
template class tbox::ListIterator< xfer::CoarsenClasses<NDIM>::Data >;
template class tbox::ListNode< xfer::CoarsenClasses<NDIM>::Data >;
}
}
