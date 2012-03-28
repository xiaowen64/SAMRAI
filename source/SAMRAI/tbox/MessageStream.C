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

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * The constructor and destructor for MessageStream.
 *
 *************************************************************************
 */

MessageStream::MessageStream(
   const size_t num_bytes,
   const StreamMode mode,
   const void *data_to_read):
   d_mode(mode),
   d_buffer(),
   d_buffer_index(0),
   d_grow_as_needed(false)
{
   TBOX_ASSERT(num_bytes >= 1);
   d_buffer.reserve(num_bytes);

   if ( mode == Read ) {
      if ( num_bytes > 0 && data_to_read == NULL ) {
         TBOX_ERROR("MessageStream::MessageStream: error:\n"
                    <<"No data_to_read was given to a Read-mode MessageStream.\n");
      }
      d_buffer.insert( d_buffer.end(),
                       static_cast<const char*>(data_to_read),
                       static_cast<const char*>(data_to_read)+num_bytes );
   }
   return;
}

MessageStream::~MessageStream()
{
}

/*
 *************************************************************************
 *
 * Print out class data if an assertion is thrown.
 *
 *************************************************************************
 */

void
MessageStream::printClassData(
   std::ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer.size() << std::endl;
   os << "Current buffer index = " << d_buffer_index << std::endl;
   os << "Pointer to buffer data = " << static_cast<const void *>(&d_buffer[0]) << std::endl;
}

}
}

#endif
