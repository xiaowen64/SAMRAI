//
// File:	OuteredgeVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 329 $
// Modified:	$Date: 2005-05-03 11:26:58 -0700 (Tue, 03 May 2005) $
// Description:	Variable class for defining outeredge centered variables
//

#ifndef included_pdat_OuteredgeVariable_C
#define included_pdat_OuteredgeVariable_C

#include "OuteredgeVariable.h"
#include "OuteredgeDataFactory.h"
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
OuteredgeVariable<DIM,TYPE>::OuteredgeVariable(
   const string &name, int depth)
:  hier::Variable<DIM>(name, new OuteredgeDataFactory<DIM,TYPE>(depth))
{
}

template <int DIM, class TYPE>
OuteredgeVariable<DIM,TYPE>::~OuteredgeVariable()
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
OuteredgeVariable<DIM,TYPE>::OuteredgeVariable(
   const OuteredgeVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template <int DIM, class TYPE>
void OuteredgeVariable<DIM,TYPE>::operator=(
   const OuteredgeVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif

