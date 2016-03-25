//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outernode/OuternodeVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
   const std::string &name, int depth)
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

