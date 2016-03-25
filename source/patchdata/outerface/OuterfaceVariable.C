//
// File:	OuterfaceVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
   const string &name, int depth)
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
