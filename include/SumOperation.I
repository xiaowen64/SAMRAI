//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
   vdst += vsrc;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
