//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/node/NodeVariable.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Variable class for defining node centered variables
//

#ifndef included_pdat_NodeVariable_C
#define included_pdat_NodeVariable_C

#include "NodeVariable.h"
#include "NodeDataFactory.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for node variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
NodeVariable<DIM,TYPE>::NodeVariable(
   const std::string &name,
   int depth, 
   bool fine_boundary_represents_var)
:  hier::Variable<DIM>(name, 
                  new NodeDataFactory<DIM,TYPE>(depth, 
                                                  // default zero ghost cells
                                                  hier::IntVector<DIM>(0),
                                                  fine_boundary_represents_var)),
   
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<int DIM, class TYPE>
NodeVariable<DIM,TYPE>::~NodeVariable()
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
NodeVariable<DIM,TYPE>::NodeVariable(
   const NodeVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void NodeVariable<DIM,TYPE>::operator=(const NodeVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
