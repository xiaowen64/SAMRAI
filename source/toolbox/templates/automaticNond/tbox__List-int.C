//
// File:	tbox__List-int.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< int > *tbox::ListNode< int >::s_free_list=0;
bool tbox::ListNode< int >::s_registered_callback=false;
#endif
template class tbox::List< int >;
template class tbox::ListIterator< int >;
template class tbox::ListNode< int >;
}
}
