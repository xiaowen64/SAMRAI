//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/Tracer.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
