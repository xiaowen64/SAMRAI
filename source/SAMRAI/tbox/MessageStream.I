/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
void *MessageStream::getBufferStart()
{
   return (void *)d_buffer;
}

SAMRAI_INLINE_KEYWORD
size_t MessageStream::getCurrentSize() const
{
   return d_current_size;
}

SAMRAI_INLINE_KEYWORD
size_t MessageStream::getCurrentIndex() const
{
   return d_buffer_index;
}

SAMRAI_INLINE_KEYWORD
void MessageStream::setCurrentIndex(
   const size_t index)
{
   d_buffer_index = index;
}

SAMRAI_INLINE_KEYWORD
void MessageStream::resetIndex()
{
   setCurrentIndex(0);
}

SAMRAI_INLINE_KEYWORD
void *
MessageStream::getPointerAndAdvanceCursor(const size_t nbytes)
{
   void* ptr = &d_buffer[d_buffer_index];
   d_buffer_index += nbytes;

   TBOX_ASSERT(d_buffer_index <= d_buffer_size);

   d_current_size = d_buffer_index;
   return ptr;
}

}
}
