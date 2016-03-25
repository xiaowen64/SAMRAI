//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/GridGeometry.I $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
