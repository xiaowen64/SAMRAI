//
// File:	IndexVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	0.1
// Revision:	$Revision: 605 $
// Modified:	$Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description:	hier::Variable class for defining irregular index variables
//

#ifndef included_pdat_IndexVariable_C
#define included_pdat_IndexVariable_C

#include "IndexVariable.h"
#include "IndexDataFactory.h"

#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for irregular index variable objects	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
IndexVariable<DIM,TYPE>::IndexVariable(const string &name)
:  hier::Variable<DIM>(name, new IndexDataFactory<DIM,TYPE>(hier::IntVector<DIM>(0))) 
                                                        // default zero ghost cells
{
}

template<int DIM, class TYPE>
IndexVariable<DIM,TYPE>::~IndexVariable()
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
IndexVariable<DIM,TYPE>::IndexVariable(
   const IndexVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   // not implemented
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void IndexVariable<DIM,TYPE>::operator=(const IndexVariable<DIM,TYPE>& foo)
{
   // not implemented
   NULL_USE(foo);
}

}
}
#endif
