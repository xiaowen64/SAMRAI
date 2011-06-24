/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Copy operation on single array data elements templated on data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

/*
 * Member functions for CopyOperation
 */

template<class TYPE>
SAMRAI_INLINE_KEYWORD
CopyOperation<TYPE>::CopyOperation()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
CopyOperation<TYPE>::~CopyOperation()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void CopyOperation<TYPE>::operator () (
   TYPE& vdst,
   const TYPE& vsrc) const
{
   vdst = vsrc;
}

}
}
