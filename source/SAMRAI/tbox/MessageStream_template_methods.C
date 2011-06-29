/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication 
 *
 ************************************************************************/

#ifndef included_tbox_MessageStream_C
#define included_tbox_MessageStream_C

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cstring>

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Pack array into message stream.
 *
 *************************************************************************
 */

template<typename DATA_TYPE>
void MessageStream::pack(
   const DATA_TYPE* data,
   unsigned int size)
{
   TBOX_ASSERT(d_mode == MessageStream::Write);

   if (data && (size > 0)) {
      const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
      void* ptr = getPointerAndAdvanceCursor(nbytes);
      memcpy(ptr, static_cast<const void *>(data), nbytes);
   }
}

/*
 *************************************************************************
 *
 * Pack single data value into message stream.
 *
 *************************************************************************
 */

template<typename DATA_TYPE>
MessageStream & MessageStream::operator << (
   const DATA_TYPE& data)
{
   TBOX_ASSERT(d_mode == MessageStream::Write);

   static const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(1);

   void* ptr = getPointerAndAdvanceCursor(nbytes);
   memcpy(ptr, static_cast<const void *>(&data), nbytes);

   return *this;
}

/*
 *************************************************************************
 *
 * Unpack array from message stream.
 *
 *************************************************************************
 */

template<typename DATA_TYPE>
void MessageStream::unpack(
   DATA_TYPE* data,
   unsigned int size)
{
   TBOX_ASSERT(d_mode == MessageStream::Read);

   if (data && (size > 0)) {
      const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
      void* ptr = getPointerAndAdvanceCursor(nbytes);
      memcpy(static_cast<void *>(data), ptr, nbytes);
   }
}

/*
 *************************************************************************
 *
 * Unpack single data value from message stream.
 *
 *************************************************************************
 */

template<typename DATA_TYPE>
MessageStream & MessageStream::operator >> (
   DATA_TYPE& data)
{
   TBOX_ASSERT(d_mode == MessageStream::Read);

   static const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(1);

   void* ptr = getPointerAndAdvanceCursor(nbytes);
   memcpy(static_cast<void *>(&data), ptr, nbytes);

   return *this;
}

}
}

#endif
