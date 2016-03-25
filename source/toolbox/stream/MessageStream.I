//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/stream/MessageStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Fixed-size message buffer used in interprocessor communication
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

#include "tbox/Utilities.h"

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

