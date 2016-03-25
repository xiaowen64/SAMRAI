//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/sundials/tbox_Pointer-CVODEModel.C $
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2002 Lawrence Livermore National Security, LLC
//

#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#include "CVODEModel.h"

#if defined(HAVE_SUNDIALS) && defined(HAVE_HYPRE)

namespace SAMRAI {

template class tbox::Pointer< CVODEModel >;

}
#endif
