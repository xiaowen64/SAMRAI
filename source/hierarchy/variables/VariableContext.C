//
// File:	VariableContext.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Simple integer id and namestring variable context
//

#ifndef included_hier_VariableContext_C
#define included_hier_VariableContext_C

#include "VariableContext.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#include "VariableContext.I"
#endif

namespace SAMRAI {
   namespace hier {

int VariableContext::s_instance_counter = 0;

/*
*************************************************************************
*                                                                       *
* The constructor copies the name of the variable context, obtains      *
* a unique instance number, and increments the number of global         *
* instances.  The destructor releases the name storage but does not     *
* decrease the instance count; instance numbers are not recycled.       *
*                                                                       *
*************************************************************************
*/

VariableContext::VariableContext(const string& name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif
   d_index = s_instance_counter++;
   d_name = name;
}

VariableContext::~VariableContext()
{
}

}
}

#endif
