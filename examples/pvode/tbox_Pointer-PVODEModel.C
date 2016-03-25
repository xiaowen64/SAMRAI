//
// File:	Pointer-PVODEModel.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2002 The Regents of the University of California
//

#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#include "PVODEModel.h"

#if defined(HAVE_PVODE) && defined(HAVE_HYPRE)

namespace SAMRAI {

template class tbox::Pointer< PVODEModel >;

}
#endif
