//
// File:	tbox__List-mblk__MultiblockRefineSchedule_NDIM___SingularityPatchX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "MultiblockRefineSchedule.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch > *tbox::ListNode< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch >::s_free_list=0;
bool tbox::ListNode< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch >::s_registered_callback=false;
#endif
template class tbox::List< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch >;
template class tbox::ListIterator< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch >;
template class tbox::ListNode< mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch >;
}
}
