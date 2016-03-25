//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/VariableContext.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
