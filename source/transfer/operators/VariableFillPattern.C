//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/xfer/VariableFillPattern.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2861 $
// Modified:	$LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description:	Abstract fill pattern class to provide interface for stencils
//

#ifndef included_xfer_VariableFillPattern_C
#define included_xfer_VariableFillPattern_C

#include "VariableFillPattern.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor                                                   *
*                                                                       *
*************************************************************************
*/
template<int DIM>
VariableFillPattern<DIM>::VariableFillPattern()
{
}

/*
*************************************************************************
*									*
* Destructor                                                            *
*									*
*************************************************************************
*/
template<int DIM>
VariableFillPattern<DIM>::~VariableFillPattern()
{
}



}
}
#endif
