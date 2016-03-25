//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/templates/special/tbox_List-tbox_TimerPointer.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
