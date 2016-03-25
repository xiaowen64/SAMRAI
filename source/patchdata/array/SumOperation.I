//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1889 $
// Modified:	$LastChangedDate: 2008-01-22 16:46:52 -0800 (Tue, 22 Jan 2008) $
// Description:	Sum operation on single array data elements templated on data type
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

/*
 * Member functions for SumOperation
 */

template <class TYPE> inline
SumOperation<TYPE>::SumOperation()
{
}

template <class TYPE> inline
SumOperation<TYPE>::~SumOperation()
{
}

template <class TYPE> inline
void SumOperation<TYPE>::operator()(TYPE& vdst, const TYPE& vsrc) const
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   vdst += vsrc;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
