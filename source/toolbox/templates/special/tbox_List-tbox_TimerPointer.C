//
// File:	tbox_List-tbox_TimerPointer.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "tbox/Timer.h"

namespace SAMRAI {
	namespace tbox {

#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< tbox::Timer* > *tbox::ListNode< tbox::Timer* >::s_free_list=0;
bool tbox::ListNode< tbox::Timer* >::s_registered_callback=false;
#endif

template class tbox::List< tbox::Timer* >;
template class tbox::ListIterator< tbox::Timer* >;
template class tbox::ListNode< tbox::Timer* >;

}
}
