//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerface/OuterfaceVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Variable class for defining outerface centered variables
//

#ifndef included_pdat_OuterfaceVariable_C
#define included_pdat_OuterfaceVariable_C

#include "OuterfaceVariable.h"
#include "OuterfaceDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for face variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OuterfaceVariable<DIM,TYPE>::OuterfaceVariable(
   const std::string &name, int depth)
:  hier::Variable<DIM>(name, new OuterfaceDataFactory<DIM,TYPE>(depth))
{
}

template<int DIM, class TYPE>
OuterfaceVariable<DIM,TYPE>::~OuterfaceVariable()
{
}

/*
*************************************************************************
*									*
* These are private and should not be used.  They are defined here	*
* because some template instantiation methods fail if some member	*
* functions are left undefined.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OuterfaceVariable<DIM,TYPE>::OuterfaceVariable(
   const OuterfaceVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void OuterfaceVariable<DIM,TYPE>::operator=(
   const OuterfaceVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
