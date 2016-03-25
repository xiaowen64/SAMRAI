//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/cell/CellVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Variable class for defining cell centered variables
//

#ifndef included_pdat_CellVariable_C
#define included_pdat_CellVariable_C

#include "CellVariable.h"
#include "CellDataFactory.h"

#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for cell variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
CellVariable<DIM,TYPE>::CellVariable(
   const std::string &name,
   int depth)
:  hier::Variable<DIM>(name, 
                  new CellDataFactory<DIM,TYPE>(depth, 
                                                  // default zero ghost cells
                                                  hier::IntVector<DIM>(0)))
{
}

template<int DIM, class TYPE>
CellVariable<DIM,TYPE>::~CellVariable()
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
CellVariable<DIM,TYPE>::CellVariable(
   const CellVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void CellVariable<DIM,TYPE>::operator=(const CellVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
