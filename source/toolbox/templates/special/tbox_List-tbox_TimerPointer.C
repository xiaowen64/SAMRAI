//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/templates/special/tbox_List-tbox_TimerPointer.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
