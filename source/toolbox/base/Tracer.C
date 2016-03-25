//
// File:        Tracer.C
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A simple call sequence tracking class
//

#include "tbox/Tracer.h"
#include "tbox/PIO.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Tracer.I"
#endif

namespace SAMRAI {
   namespace tbox {

ostream* Tracer::s_stream = &plog;

}
}
