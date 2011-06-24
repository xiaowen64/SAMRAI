/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   A simple call sequence tracking class 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Tracer.h"
#include "SAMRAI/tbox/PIO.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/Tracer.I"
#endif

namespace SAMRAI {
namespace tbox {

std::ostream * Tracer::s_stream = &plog;

}
}
