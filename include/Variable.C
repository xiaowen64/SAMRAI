//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/Variable.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
   const std::string &name,
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
