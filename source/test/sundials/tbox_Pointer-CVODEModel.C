/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:
 *
 ************************************************************************/
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Pointer.C"
#include "CVODEModel.h"

#if defined(HAVE_SUNDIALS) && defined(HAVE_HYPRE)

namespace SAMRAI {

template class tbox::Pointer<CVODEModel>;

}
#endif
