//
// File:	Variable.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Base class for application-level variables
//

#ifndef included_hier_Variable_C
#define included_hier_Variable_C

#include "Variable.h"

#ifdef DEBUG_NO_INLINE
#include "Variable.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> int Variable<DIM>::s_instance_counter = 0;

/*
*************************************************************************
*									*
* The constructor copies the name of the variable, obtains a unique	*
* instance number, and increments the number of global instances.	*
* The destructor releases the name storage but does not decrease the	*
* instance count, since instance numbers are never recycled.		*
*									*
*************************************************************************
*/

template<int DIM>  Variable<DIM>::Variable(
   const string &name,
   const tbox::Pointer< PatchDataFactory<DIM> > factory)
{
   d_instance = s_instance_counter++;
   d_name = name;
   d_factory = factory;
}

template<int DIM>  Variable<DIM>::~Variable()
{
}

}
}
#endif
