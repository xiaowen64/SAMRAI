//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/Tracer.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: A simple call sequence tracking class
//

#include "tbox/Tracer.h"
#include "tbox/PIO.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Tracer.I"
#endif

namespace SAMRAI {
   namespace tbox {

std::ostream* Tracer::s_stream = &plog;

}
}
