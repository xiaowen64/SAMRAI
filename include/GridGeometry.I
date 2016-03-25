//
// File:	GridGeometry.I
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
