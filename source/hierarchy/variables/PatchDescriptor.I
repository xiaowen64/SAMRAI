//
// File:	PatchDescriptor.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 407 $
// Modified:	$Date: 2005-06-01 09:49:54 -0700 (Wed, 01 Jun 2005) $
// Description:	Factory class for patch data objects that live on a patch
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
tbox::Pointer< PatchDataFactory<DIM> >
PatchDescriptor<DIM>::getPatchDataFactory(const int id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (id >= 0) && (id < d_max_number_registered_components) );
#endif
   return(d_factories[id]);
}

template<int DIM> inline
int PatchDescriptor<DIM>::getMaxNumberRegisteredComponents() const
{
   return(d_max_number_registered_components);
}

template<int DIM> inline
const string& PatchDescriptor<DIM>::mapIndexToName(const int id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (id >= 0) && (id < d_max_number_registered_components) );
#endif
   return(d_names[id]);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
