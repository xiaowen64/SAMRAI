//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/stream/IOStream.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Wrapper header file for standard IO stream classes
//

#ifndef included_Stream
#define included_Stream

#include "SAMRAI_config.h"

#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_iomanip
#define included_iomanip
#include <iomanip>
#endif

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#endif

#ifndef LACKS_STRSTREAM
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif


#endif
