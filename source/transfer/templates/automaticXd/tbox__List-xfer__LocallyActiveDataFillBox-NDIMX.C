//
// File:	tbox__List-xfer__LocallyActiveDataFillBox-NDIMX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "LocallyActiveDataFillBox.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< xfer::LocallyActiveDataFillBox< NDIM > > *tbox::ListNode< xfer::LocallyActiveDataFillBox< NDIM > >::s_free_list=0;
bool tbox::ListNode< xfer::LocallyActiveDataFillBox< NDIM > >::s_registered_callback=false;
#endif
template class tbox::List< xfer::LocallyActiveDataFillBox< NDIM > >;
template class tbox::ListIterator< xfer::LocallyActiveDataFillBox< NDIM > >;
template class tbox::ListNode< xfer::LocallyActiveDataFillBox< NDIM > >;
}
}
