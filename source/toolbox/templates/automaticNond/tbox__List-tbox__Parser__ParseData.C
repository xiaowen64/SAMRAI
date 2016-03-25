//
// File:	tbox__List-tbox__Parser__ParseData.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 1.32 
// Modified:	$Date: 2003/01/22 01:12:43 
// Description:	Automatically generated template file
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "tbox/Parser.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< tbox::Parser::ParseData > *tbox::ListNode< tbox::Parser::ParseData >::s_free_list=0;
bool tbox::ListNode< tbox::Parser::ParseData >::s_registered_callback=false;
#endif
template class tbox::List< tbox::Parser::ParseData >;
template class tbox::ListIterator< tbox::Parser::ParseData >;
template class tbox::ListNode< tbox::Parser::ParseData >;
}
}
