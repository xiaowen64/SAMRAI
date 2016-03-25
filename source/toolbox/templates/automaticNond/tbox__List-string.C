//
// File:	tbox__List-string.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include <string>
using namespace std;

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< string > *tbox::ListNode< string >::s_free_list=0;
bool tbox::ListNode< string >::s_registered_callback=false;
#endif
template class tbox::List< string >;
template class tbox::ListIterator< string >;
template class tbox::ListNode< string >;
}
}
