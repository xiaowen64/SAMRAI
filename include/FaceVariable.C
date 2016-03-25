//
// File:	FaceVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining face centered variables
//

#ifndef included_pdat_FaceVariable_C
#define included_pdat_FaceVariable_C

#include "FaceVariable.h"
#include "FaceDataFactory.h"
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
FaceVariable<DIM,TYPE>::FaceVariable(
   const string &name, 
   int depth, 
   const bool fine_boundary_represents_var)
:  hier::Variable<DIM>(name, 
                  new FaceDataFactory<DIM,TYPE>(depth, 
                                                  // default zero ghost cells
                                                  hier::IntVector<DIM>(0),
                                                  fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<int DIM, class TYPE>
FaceVariable<DIM,TYPE>::~FaceVariable()
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
FaceVariable<DIM,TYPE>::FaceVariable(
   const FaceVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void FaceVariable<DIM,TYPE>::operator=(const FaceVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
