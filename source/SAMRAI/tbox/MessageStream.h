/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/

#ifndef included_tbox_MessageStream
#define included_tbox_MessageStream

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cstring>
#include <iostream>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class to provide buffers for communication of data.
 *
 * MessageStream provides a fixed size message buffer that can
 * hold data of any type.  It is used by communication routines in the
 * Schedule class.
 *
 * @see tbox::Schedule
 */

class MessageStream
{
public:
   /*!
    * @brief Enumeration to identify if a buffer is being used to read or
    * write data.
    */
   enum StreamMode { Read, Write };

   /*!
    * @brief Create a message stream of the specified size and mode
    *
    * @param[in] bytes   Number of bytes in the stream.
    * @param[in] mode    MessageStream::Read or MessageStream::Write.
    */
   MessageStream(
      const size_t bytes,
      const StreamMode mode);

   /*!
    * Destructor for a message stream.
    */
   ~MessageStream();

   /*!
    * @brief Static method to get amount of message stream space needed to
    * communicate data type indicated by template parameter.
    *
    * IMPORTANT:  All size information given to the message stream should
    * be based on values returned by this method.
    *
    * TODO:  Implementation should be moved out of header?  If we do this,
    * Then we need to create another implementation file to include in this
    * header.  I don't think it's worth it. RDH
    *
    * @return The number of bytes for num_items of type DATA_TYPE.
    *
    * @param[in] num_items
    */
   template<typename DATA_TYPE>
   static unsigned int getSizeof(
      unsigned int num_items = 1)
   {
      return num_items * static_cast<unsigned int>(sizeof(DATA_TYPE));
   }

   /*!
    * @brief Return a pointer to the start of the message buffer.
    */
   void *
   getBufferStart();

   /*!
    * @brief Return the current size of the buffer in bytes.
    */
   size_t
   getCurrentSize() const;

   /*!
    * @brief Return the current index into the buffer.
    */
   size_t
   getCurrentIndex() const;

   /*!
    * @brief Set the current index into the buffer.  Further packing/unpacking
    * will begin at this new location.
    */
   void
   setCurrentIndex(
      const size_t index);

   /*!
    * @brief Reset the index to the beginning of the buffer.  This is the same
    * as setting the buffer index to zero via setCurrentIndex().
    */
   void
   resetIndex();

   /*!
    * @brief Pack a single data item into message stream.
    *
    * @param[in] data  Single item of type DATA_TYPE to be copied
    * into the stream.
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator << (
      const DATA_TYPE& data)
   {
      TBOX_ASSERT(d_mode == MessageStream::Write);

      static const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(1);

      void* ptr = getPointerAndAdvanceCursor(nbytes);
      memcpy(ptr, static_cast<const void *>(&data), nbytes);

      return *this;
   }

   /*!
    * @brief Pack an array of data items into message stream.
    *
    * @param[in] data  Pointer to an array of data of type DATA_TYPE
    *                  to be copied into the stream.
    * @param[in] size  Number of items to pack.
    */
   template<typename DATA_TYPE>
   void
   pack(
      const DATA_TYPE* data,
      unsigned int size = 1)
   {
      TBOX_ASSERT(d_mode == MessageStream::Write);

      if (data && (size > 0)) {
         const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         void* ptr = getPointerAndAdvanceCursor(nbytes);
         memcpy(ptr, static_cast<const void *>(data), nbytes);
      }
   }

   /*!
    * @brief Unpack a single data item from message stream.
    *
    * @param[out] data  Single item of type DATA_TYPE that will be
    *                   copied from the stream.
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator >> (
      DATA_TYPE& data)
   {
      TBOX_ASSERT(d_mode == MessageStream::Read);

      static const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(1);

      void* ptr = getPointerAndAdvanceCursor(nbytes);
      memcpy(static_cast<void *>(&data), ptr, nbytes);

      return *this;
   }

   /*!
    * @brief Unpack an array of data items from message stream.
    *
    * @param[out] data  Pointer to an array of data of type DATA_TYPE
    *                   that will receive data copied from
    *                   the stream.
    * @param[out] size  Number of items that will be copied.
    */
   template<typename DATA_TYPE>
   void
   unpack(
      DATA_TYPE * data,
      unsigned int size = 1)
   {
      TBOX_ASSERT(d_mode == MessageStream::Read);

      if (data && (size > 0)) {
         const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         void* ptr = getPointerAndAdvanceCursor(nbytes);
         memcpy(static_cast<void *>(data), ptr, nbytes);
      }
   }

   /*!
    * @brief Print out internal object data.
    *
    * @param[out] os  Output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

private:
   /*!
    * @brief  Helper function to get pointer into buffer at current
    * position and advance buffer position by given number of bytes.
    *
    * @param[in]  nbytes
    */
   void *
   getPointerAndAdvanceCursor(
      const size_t nbytes);

   MessageStream(
      const MessageStream&);            // not implemented
   void
   operator = (
      const MessageStream&);            // not implemented

   /*!
    * @brief  Read/write mode of the stream.
    */
   const StreamMode d_mode;

   /*!
    * Number of bytes allocated in the buffer.
    */
   size_t d_buffer_size;

   /*!
    * Number of bytes currently being used in the buffer.
    */
   size_t d_current_size;

   /*!
    * Current index into the buffer used when traversing.
    */
   size_t d_buffer_index;

   /*!
    * The buffer for the streamed data.
    */
   char* d_buffer;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/MessageStream.I"
#endif

#endif
