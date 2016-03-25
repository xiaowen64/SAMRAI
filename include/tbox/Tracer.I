//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/Tracer.I $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
