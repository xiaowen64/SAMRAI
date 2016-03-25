//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/edge/EdgeVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Variable class for defining edge centered variables
//

#ifndef included_pdat_EdgeVariable_C
#define included_pdat_EdgeVariable_C

#include "EdgeVariable.h"
#include "EdgeDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for edge variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
EdgeVariable<DIM,TYPE>::EdgeVariable(
   const std::string &name,
   int depth, 
   const bool fine_boundary_represents_var)
:  hier::Variable<DIM>(name, 
                  new EdgeDataFactory<DIM,TYPE>(depth, 
                                                  // default zero ghost cells
                                                  hier::IntVector<DIM>(0), 
                                                  fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var) 
{
}

template<int DIM, class TYPE>
EdgeVariable<DIM,TYPE>::~EdgeVariable()
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
EdgeVariable<DIM,TYPE>::EdgeVariable(
   const EdgeVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void EdgeVariable<DIM,TYPE>::operator=(const EdgeVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
