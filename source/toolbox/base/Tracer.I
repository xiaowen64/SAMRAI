//
// File:        Tracer.I
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A simple call sequence tracking class
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
Tracer::Tracer(const string& message)
{
   d_message = message;
   if (s_stream) {
      (*s_stream) << "Entering " << d_message << endl << flush;
   }
}

inline
Tracer::~Tracer()
{
   if (s_stream) {
      (*s_stream) << "Exiting " << d_message << endl << flush;
   }
}

inline
void Tracer::setTraceStream(ostream* stream)
{
   s_stream = stream;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
