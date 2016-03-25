//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/stream/XDRStream.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Stream class that converts into XDR for portable communication
//

#include "tbox/XDRStream.h"
#include <string>
using namespace std;
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/XDRStream.I"
#endif

namespace SAMRAI {
   namespace tbox {



/*
*************************************************************************
*									*
* Some macros to simplify the implementation for the XDR streams	*
*									*
*************************************************************************
*/

#ifdef HAVE_XDR

#define XDR_PACK_OPAQUE(m_data,m_size)					\
   ((d_xdr_stream->x_op != XDR_ENCODE) ||				\
    (!xdr_opaque(d_xdr_stream, (caddr_t) m_data, m_size)))
#define XDR_UNPACK_OPAQUE(m_data,m_size)				\
   ((d_xdr_stream->x_op != XDR_DECODE) ||				\
    (!xdr_opaque(d_xdr_stream, (caddr_t) m_data, m_size)))
#define XDR_PACK_VECTOR(m_data,m_size,m_type)				\
   ((d_xdr_stream->x_op != XDR_ENCODE) ||				\
    (!xdr_vector(d_xdr_stream, (char *) m_data, m_size,			\
                 sizeof(m_type), (xdrproc_t) xdr_##m_type)))
#define XDR_UNPACK_VECTOR(m_data,m_size,m_type)				\
   ((d_xdr_stream->x_op != XDR_DECODE) ||				\
    (!xdr_vector(d_xdr_stream, (char *) m_data, m_size,			\
                 sizeof(m_type), (xdrproc_t) xdr_##m_type)))

#else

#define XDR_PACK_OPAQUE(m_data,m_size) 0
#define XDR_UNPACK_OPAQUE(m_data,m_size) 0
#define XDR_PACK_VECTOR(m_data,m_size,m_type) 0
#define XDR_UNPACK_VECTOR(m_data,m_size,m_type) 0

#endif


/*
*************************************************************************
*									*
* The virtual destructor for XDRStream does nothing.		*
*									*
*************************************************************************
*/

XDRStream::~XDRStream()
{
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* into character arrays and then packed using the character routines.	*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const bool& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(bool& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const bool *data, const int n)
{
   char *flags = new char[n];
   for (int i = 0; i < n; i++) {
      flags[i] = (data[i] ? 1 : 0);
   }
   if (XDR_PACK_OPAQUE(flags, n)) {
      TBOX_ERROR("XDRStream: Error in encoding bool...\n");
   }
   delete [] flags;
}

void XDRStream::unpack(bool *data, const int n)
{
   char *flags = new char[n];
   if (XDR_UNPACK_OPAQUE(flags, n)) {
      TBOX_ERROR("XDRStream: Error in decoding bool...\n");
   }
   for (int i = 0; i < n; i++) {
      data[i] = (flags[i] ? true : false);
   }
   delete [] flags;
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const char& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(char& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const char *data, const int n)
{
   if (XDR_PACK_OPAQUE(data, n)) {
      TBOX_ERROR("XDRStream: Error in encoding char...\n");
   }
}

void XDRStream::unpack(char *data, const int n)
{
   if (XDR_UNPACK_OPAQUE(data, n)) {
      TBOX_ERROR("XDRStream: Error in decoding char...\n");
   }
}

void XDRStream::writeString(const char *data)
{
#ifdef HAVE_XDR
   if (!xdr_string(d_xdr_stream, (char **) &data, strlen(data))) {
      TBOX_ERROR("XDRStream: Error in writing string...\n");
   }
#endif

}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex		*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const dcomplex& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(dcomplex& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(sizeof(dcomplex) == 2*sizeof(double));
#endif
   if (XDR_PACK_VECTOR((double *) data, 2*n, double)) {
      TBOX_ERROR("XDRStream: Error in encoding double complex...\n");
   }
}

void XDRStream::unpack(dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(sizeof(dcomplex) == 2*sizeof(double));
#endif
   if (XDR_UNPACK_VECTOR((double *) data, 2*n, double)) {
      TBOX_ERROR("XDRStream: Error in decoding double complex...\n");
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const double& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(double& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const double *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, double)) {
      TBOX_ERROR("XDRStream: Error in encoding double...\n");
   }
}

void XDRStream::unpack(double *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, double)) {
      TBOX_ERROR("XDRStream: Error in decoding double...\n");
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const float& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(float& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const float *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, float)) {
      TBOX_ERROR("XDRStream: Error in encoding float...\n");
   }
}

void XDRStream::unpack(float *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, float)) {
      TBOX_ERROR("XDRStream: Error in decoding float...\n");
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const int& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(int& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const int *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, int)) {
      TBOX_ERROR("XDRStream: Error in encoding integer...\n");
   }
}

void XDRStream::unpack(int *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, int)) {
      TBOX_ERROR("XDRStream: Error in decoding integer...\n");
   }
}

}
}
