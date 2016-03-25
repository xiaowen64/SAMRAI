//
// File:	Pointer-MultiblockTester.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2005 The Regents of the University of California
//

#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#include "MultiblockTester.h"

using namespace SAMRAI;

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION
template class tbox::Pointer< MultiblockTester >;
#endif

