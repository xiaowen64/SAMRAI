//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/Tracer.h $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1818 $
// Modified:    $LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
// Description: A simple call sequence tracking class
//

#ifndef included_tbox_Tracer
#define included_tbox_Tracer

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#include <iostream>
#define included_iostream
#endif

#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class Tracer allows one to trace entrances and exits of
 * class member functions. An example usage is:
 * \verbatim
 * #include "tbox/Tracer.h"
 * ....
 * void MyClass::myClassMemberFunction() 
 * {
 *    Tracer t("MyClass::myClassMemberFunction");
 *     ....
 * }
 * \endverbatim
 * When the function `myClassMemberFunction' is called, a tracer
 * object  local to the function scope is created and the message
 * {\em Entering MyClass::myClassMemberFunction} is sent to the
 * tracer output stream. Upon exiting the function, the tracer object
 * is destroyed and {\em Exiting MyClass::myClassMemberFunction}
 * is sent to the tracer output stream.  By default, the tracer
 * class sends data to the parallel log stream, although the default
 * output stream can be changed through a call to static member function
 * setTraceStream(), which will set the tracer output stream for all
 * subsequent calls.
 */

class Tracer 
{
public:
   /**
    * The constructor for Tracer prints ``Entering \<message\>''
    * to the tracer output stream.
    */
   Tracer(const std::string& message);

   /**
    * The destructor for Tracer prints ``Exiting \<message\>''
    * to the tracer output stream.
    */
    ~Tracer();

   /**
    * Set the tracer output stream for all tracer output.  By default,
    * this is set to the parallel log stream plog.  If this argument is
    * NULL, then all output to trace streams is disabled.
    */
   static void setTraceStream(std::ostream* stream);

private:
   std::string d_message;
   static std::ostream* s_stream;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Tracer.I"
#endif
#endif
