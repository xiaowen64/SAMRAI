//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/Tracer.I $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: A simple call sequence tracking class
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
Tracer::Tracer(const std::string& message)
{
   d_message = message;
   if (s_stream) {
      (*s_stream) << "Entering " << d_message << std::endl << std::flush;
   }
}

inline
Tracer::~Tracer()
{
   if (s_stream) {
      (*s_stream) << "Exiting " << d_message << std::endl << std::flush;
   }
}

inline
void Tracer::setTraceStream(std::ostream* stream)
{
   s_stream = stream;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
