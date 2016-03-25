//
// File:	IOStream.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Wrapper header file for standard IO stream classes
//

#ifndef included_Stream
#define included_Stream

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

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
using namespace std;


#endif
