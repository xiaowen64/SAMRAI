/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Sum operation on single array data elements templated on data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

/*
 * Member functions for SumOperation
 */

template<class TYPE>
SAMRAI_INLINE_KEYWORD
SumOperation<TYPE>::SumOperation()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
SumOperation<TYPE>::~SumOperation()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void SumOperation<TYPE>::operator () (
   TYPE& vdst,
   const TYPE& vsrc) const
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   vdst += vsrc;
}

}
}
