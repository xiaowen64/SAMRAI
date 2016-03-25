//
// File:	tbox__List-appu__VisItDataWriter_NDIM___VisItItemX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "VisItDataWriter.h"

namespace SAMRAI {
   namespace tbox {
#ifdef HAVE_HDF5
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< appu::VisItDataWriter<NDIM>::VisItItem > *tbox::ListNode< appu::VisItDataWriter<NDIM>::VisItItem >::s_free_list=0;
bool tbox::ListNode< appu::VisItDataWriter<NDIM>::VisItItem >::s_registered_callback=false;
#endif
template class tbox::List< appu::VisItDataWriter<NDIM>::VisItItem >;
template class tbox::ListIterator< appu::VisItDataWriter<NDIM>::VisItItem >;
template class tbox::ListNode< appu::VisItDataWriter<NDIM>::VisItItem >;
#endif
}
}
