//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/GridGeometry.I $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Base class for geometry management in AMR hierarchy
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
const BoxArray<DIM>& GridGeometry<DIM>::getPhysicalDomain() const
{
   return (d_physical_domain);
}

template<int DIM> inline
bool GridGeometry<DIM>::getDomainIsSingleBox() const
{
   return (d_domain_is_single_box);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
