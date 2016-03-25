//
// File:	tbox__List-hier__Box-NDIMX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "Box.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< hier::Box< NDIM > > *tbox::ListNode< hier::Box< NDIM > >::s_free_list=0;
bool tbox::ListNode< hier::Box< NDIM > >::s_registered_callback=false;
#endif
template class tbox::List< hier::Box< NDIM > >;
template class tbox::ListIterator< hier::Box< NDIM > >;
template class tbox::ListNode< hier::Box< NDIM > >;
}
}
