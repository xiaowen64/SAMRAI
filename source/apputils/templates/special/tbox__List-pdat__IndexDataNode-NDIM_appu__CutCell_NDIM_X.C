//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/templates/special/tbox__List-pdat__IndexDataNode-NDIM_appu__CutCell_NDIM_X.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2005/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "IndexData.h"
#include "CutCell.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > > *tbox::ListNode< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > >::s_free_list=0;
bool tbox::ListNode< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > >::s_registered_callback=false;
#endif
template class tbox::List< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > >;
template class tbox::ListIterator< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > >;
template class tbox::ListNode< pdat::IndexDataNode< NDIM,appu::CutCell<NDIM> > >;
}
}
