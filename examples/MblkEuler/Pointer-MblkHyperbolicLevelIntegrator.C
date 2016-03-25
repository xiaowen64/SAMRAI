//
// File:	Pointer-MultiblockTester.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
//

#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#include "MblkHyperbolicLevelIntegrator.h"

using namespace SAMRAI;

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION
template class tbox::Pointer< MblkHyperbolicLevelIntegrator >;
#endif

