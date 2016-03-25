//
// File:	tbox__List-xfer__RefineClasses_NDIM___DataX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "RefineClasses.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< xfer::RefineClasses<NDIM>::Data > *tbox::ListNode< xfer::RefineClasses<NDIM>::Data >::s_free_list=0;
bool tbox::ListNode< xfer::RefineClasses<NDIM>::Data >::s_registered_callback=false;
#endif
template class tbox::List< xfer::RefineClasses<NDIM>::Data >;
template class tbox::ListIterator< xfer::RefineClasses<NDIM>::Data >;
template class tbox::ListNode< xfer::RefineClasses<NDIM>::Data >;
}
}
