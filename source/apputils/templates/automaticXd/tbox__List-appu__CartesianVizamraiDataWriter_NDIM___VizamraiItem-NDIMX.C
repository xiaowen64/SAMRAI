//
// File:	tbox__List-appu__CartesianVizamraiDataWriter_NDIM___VizamraiItem-NDIMX.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "CartesianVizamraiDataWriter.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > > *tbox::ListNode< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > >::s_free_list=0;
bool tbox::ListNode< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > >::s_registered_callback=false;
#endif
template class tbox::List< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > >;
template class tbox::ListIterator< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > >;
template class tbox::ListNode< appu::CartesianVizamraiDataWriter<NDIM>::VizamraiItem< NDIM > >;
}
}
