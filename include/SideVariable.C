//
// File:	SideVariable.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining side centered variables
//

#ifndef included_pdat_SideVariable_C
#define included_pdat_SideVariable_C

#include "SideVariable.h"
#include "SideDataFactory.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#define ALL_DIRECTIONS (-1)

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Constructor and destructor for side variable objects			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
SideVariable<DIM,TYPE>::SideVariable(
   const string &name, 
   int depth, 
   bool fine_boundary_represents_var,
   int direction)
:  hier::Variable<DIM>(name, 
                  new SideDataFactory<DIM,TYPE>(depth, 
                                                  // default zero ghost cells
                                                  hier::IntVector<DIM>(0), 
                                                  fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((direction >= ALL_DIRECTIONS) && (direction < DIM));
#endif
   d_directions = hier::IntVector<DIM>(1);
   if ( (direction != ALL_DIRECTIONS) ) {
      for (int id = 0; id < DIM; id++) {
         d_directions(id) = ( (direction == id) ? 1 : 0 );  
      }  
      setPatchDataFactory(new SideDataFactory<DIM,TYPE>(depth,
                                                          hier::IntVector<DIM>(0),
                                                          d_directions,
                                                          fine_boundary_represents_var));
   }
}

template<int DIM, class TYPE>
SideVariable<DIM,TYPE>::~SideVariable()
{
}

template<int DIM, class TYPE>
const hier::IntVector<DIM>& SideVariable<DIM,TYPE>::getDirectionVector() const
{
   return(d_directions);
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
SideVariable<DIM,TYPE>::SideVariable(
   const SideVariable<DIM,TYPE>& foo)
:  hier::Variable<DIM>(NULL, NULL)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void SideVariable<DIM,TYPE>::operator=(const SideVariable<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
