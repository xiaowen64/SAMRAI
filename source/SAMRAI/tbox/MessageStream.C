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

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/MessageStream.I"
#endif

#include <cstring>

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
   const size_t bytes,
   const StreamMode mode):
   d_mode(mode),
   d_buffer_size(bytes),
   d_current_size(0),
   d_buffer_index(0)
{
   TBOX_ASSERT(d_buffer_size >= 1);

   d_buffer = new char[d_buffer_size];
}

MessageStream::~MessageStream()
{
   delete[] d_buffer;
}

/*
 *************************************************************************
 *
 * Print out class data if an assertion is thrown.
 *
 *************************************************************************
 */

void MessageStream::printClassData(
   std::ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer_size << std::endl;
   os << "Current buffer size = " << d_current_size << std::endl;
   os << "Current buffer index = " << d_buffer_index << std::endl;
   os << "Pointer to buffer data = " << (void *)d_buffer << std::endl;
}

}
}

#endif
