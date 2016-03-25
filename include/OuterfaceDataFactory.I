//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerface/OuterfaceDataFactory.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Factory class for creating outerface data objects
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int OuterfaceDataFactory<DIM,TYPE>::getDefaultDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
void OuterfaceDataFactory<DIM,TYPE>::setDefaultDepth(const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   d_depth = depth;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
