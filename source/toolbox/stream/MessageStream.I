//
// File:	MessageStream.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Fixed-size message buffer used in interprocessor communication
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
   namespace tbox {


inline
void MessageStream::useXDR(const bool flag)
{
#ifdef HAVE_XDR
   s_use_xdr_translation = flag;
#else
   if ( flag == true )
      TBOX_ERROR("MessageStream: Error in XDR was not used in this compile...\n");
#endif
}

inline
void *MessageStream::getBufferStart()
{
   return((void *) d_buffer);
}

inline
int MessageStream::getCurrentSize() const
{
   return(d_current_size);
}

inline
int MessageStream::getCurrentIndex() const
{
   return(d_buffer_index);
}

inline
void MessageStream::setCurrentIndex(const int index)
{
   d_buffer_index = index;
	
   // SGS this is broken, in fact this whole class 
   // looks broken, if xdr is not used how is d_buffer_index 
   // moved?	

#ifdef HAVE_XDR
   if (d_use_xdr) {
#ifndef LACKS_PROPER_XDR_HEADER
      xdr_setpos(&d_xdr_stream, d_buffer_index);
#else
      (*(void(*)(XDR*,unsigned int))(d_xdr_stream.x_ops->x_setpostn))(
         &d_xdr_stream, d_buffer_index);
#endif
   }
#endif

}

inline
void MessageStream::resetIndex()
{
   setCurrentIndex(0);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

