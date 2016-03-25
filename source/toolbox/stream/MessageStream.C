//
// File:	MessageStream.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Fixed-size message buffer used in interprocessor communication
//

#include "tbox/MessageStream.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "tbox/MessageStream.I"
#endif

namespace SAMRAI {
   namespace tbox {


bool MessageStream::s_use_xdr_translation = false;

/*
*************************************************************************
*									*
* The constructor and destructor for MessageStream.  The first	*
* form of the constructor uses the class-wide flag to determine whether	*
* XDR will be used.  The second constructor uses the use_xdr boolean	*
* argument.  If XDR translation is to be used, then create an XDR	*
* stream object	and bind it to the XDR manager object.			*
*									*
*************************************************************************
*/

MessageStream::MessageStream(const int bytes, const StreamMode mode)
{
   d_buffer_size  = bytes;
   d_current_size = 0;
   d_buffer_index = 0;
   d_use_xdr      = s_use_xdr_translation;
   d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
   if (d_use_xdr) {
      xdr_op xop = ((mode==MessageStream::Read) ? XDR_DECODE : XDR_ENCODE);
      xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
      d_xdr_manager.setXDRStream(&d_xdr_stream);
   }
#endif

}

MessageStream::MessageStream(const int bytes,
                                      const StreamMode mode,
                                      const bool use_xdr)
{
   d_buffer_size  = bytes;
   d_current_size = 0;
   d_buffer_index = 0;
   d_use_xdr      = use_xdr;
   d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
   if (d_use_xdr) {
      xdr_op xop = ((mode==MessageStream::Read) ? XDR_DECODE : XDR_ENCODE);
      xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
      d_xdr_manager.setXDRStream(&d_xdr_stream);
   }
#endif

}

MessageStream::~MessageStream()
{
#ifdef HAVE_XDR
   if (d_use_xdr) {
#ifndef LACKS_PROPER_XDR_HEADER
      xdr_destroy(&d_xdr_stream);
#else
      if (d_xdr_stream.x_ops->x_destroy)
         (*(void(*)(XDR*))(d_xdr_stream.x_ops->x_destroy))(&d_xdr_stream);
#endif
   }
#endif

   delete [] d_buffer;
}

/*
*************************************************************************
*									*
* Print out class data if an exception is thrown.			*
*									*
*************************************************************************
*/

void MessageStream::printClassData(ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer_size << endl;
   os << "Current buffer size = " << d_current_size << endl;
   os << "Current buffer index = " << d_buffer_index << endl;
   os << "Pointer to buffer data = " << (void *) d_buffer << endl;
   os << "Using XDR translation = " << (d_use_xdr ? "true" : "false") << endl;
}

/*
*************************************************************************
*									*
* Packing/unpacking helper functions and macros.  The member function	*
* getPointerAndAdvanceCursor() returns a pointer to buffer space and	*
* advances internal pointers to reflect the allocated buffers space.	*
* The two macros given below simplify packing and unpacking for the	*
* numerous member functions below.					*
*									*
*************************************************************************
*/

void *MessageStream::getPointerAndAdvanceCursor(const int bytes)
{
   void *ptr = &d_buffer[d_buffer_index];
   d_buffer_index += bytes;
   if (d_buffer_index > d_current_size) {
      d_current_size = d_buffer_index;
      if (d_buffer_index > d_buffer_size) {
         TBOX_ERROR("MessageStream: Stream overrun of buffer...\n");
      }
   }
   return(ptr);
}

/*
*************************************************************************
*									*
* The following macros are used by all of the standard data types	*]
* except bool for packing and unpacking the data buffer.		*
*									*
*************************************************************************
*/

#ifdef HAVE_XDR

#define PACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   if (d_use_xdr) {							\
      d_xdr_manager.pack(m_data, m_size);				\
   } else {								\
      memcpy(ptr, (void *) m_data, m_bytes);				\
   }                                                                    \
   } while (0)

#define UNPACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   if (d_use_xdr) {							\
      d_xdr_manager.unpack(m_data, m_size);				\
   } else {								\
      memcpy((void *) m_data, ptr, m_bytes);				\
   }                                                                    \
   } while (0)

#else

#define PACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   memcpy(ptr, (void *) m_data, m_bytes);                               \
   } while (0)
   
#define UNPACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   memcpy((void *) m_data, ptr, m_bytes);                               \
   } while (0)

#endif


/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* either using XDR or by converting into character arrays.		*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const bool& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(bool& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const bool *data, const int n)
{
   const int bytes = AbstractStream::sizeofBool(n);
   void *ptr = getPointerAndAdvanceCursor(bytes);
   if (d_use_xdr) {
#ifdef HAVE_XDR
      d_xdr_manager.pack(data, n);
#endif
   } else {
      char *c_ptr = (char *) ptr;
      for (int i = 0; i < n; i++) {
         c_ptr[i] = (data[i] ? 1 : 0);
      }
   }
}

void MessageStream::unpack(bool *data, const int n)
{
   const int bytes = AbstractStream::sizeofBool(n);
   void *ptr = getPointerAndAdvanceCursor(bytes);
   if (d_use_xdr) {
#ifdef HAVE_XDR
      d_xdr_manager.unpack(data, 1);
#endif
   } else {
      const char *c_ptr = (const char *) ptr;
      for (int i = 0; i < n; i++) {
         data[i] = (c_ptr[i] ? true : false);
      }
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const char& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(char& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const char *data, const int n)
{
   const int bytes = AbstractStream::sizeofChar(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(char *data, const int n)
{
   const int bytes = AbstractStream::sizeofChar(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex		*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const dcomplex& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(dcomplex& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   const int bytes = AbstractStream::sizeofDoubleComplex(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   const int bytes = AbstractStream::sizeofDoubleComplex(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const double& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(double& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const double *data, const int n)
{
   const int bytes = AbstractStream::sizeofDouble(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(double *data, const int n)
{
   const int bytes = AbstractStream::sizeofDouble(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const float& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(float& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const float *data, const int n)
{
   const int bytes = AbstractStream::sizeofFloat(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(float *data, const int n)
{
   const int bytes = AbstractStream::sizeofFloat(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const int& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(int& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const int *data, const int n)
{
   const int bytes = AbstractStream::sizeofInt(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(int *data, const int n)
{
   const int bytes = AbstractStream::sizeofInt(n);
   UNPACK(data, n, bytes);
}


}
}
