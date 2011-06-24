/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   A simple call sequence tracking class 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
Tracer::Tracer(
   const std::string& message)
{
   d_message = message;
   if (s_stream) {
      (*s_stream) << "Entering " << d_message << std::endl << std::flush;
   }
}

SAMRAI_INLINE_KEYWORD
Tracer::~Tracer()
{
   if (s_stream) {
      (*s_stream) << "Exiting " << d_message << std::endl << std::flush;
   }
}

SAMRAI_INLINE_KEYWORD
void Tracer::setTraceStream(
   std::ostream* stream)
{
   s_stream = stream;
}

}
}
