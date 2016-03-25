//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/VariableContext.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Simple integer id and namestring variable context
//

#ifndef included_hier_VariableContext_C
#define included_hier_VariableContext_C

#include "VariableContext.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
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

VariableContext::VariableContext(const std::string& name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
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
