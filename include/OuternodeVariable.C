//
// File:	OuternodeVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Variable<DIM> class for defining outernode centered variables
//

#ifndef included_pdat_OuternodeVariable_C
#define included_pdat_OuternodeVariable_C

#include "OuternodeVariable.h"
#include "OuternodeDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for side variable objects                  *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuternodeVariable<DIM,TYPE>::OuternodeVariable(
   const string &name, int depth)
:  hier::Variable<DIM>(name, new OuternodeDataFactory<DIM,TYPE>(depth))
{
}

template <int DIM, class TYPE>
OuternodeVariable<DIM,TYPE>::~OuternodeVariable()
{
}

/*
*************************************************************************
*                                                                       *
* These are private and should not be used.  They are defined here      *
* because some template instantiation methods fail if some member       *
* functions are left undefined.                                         *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuternodeVariable<DIM,TYPE>::OuternodeVariable(
   const OuternodeVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template <int DIM, class TYPE>
void OuternodeVariable<DIM,TYPE>::operator=(
   const OuternodeVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif

