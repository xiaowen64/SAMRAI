//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchGeometry.I $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2147 $
// Modified:	$LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: Base class for geometry management on patches
//



#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >* PatchGeometry<DIM>::getPatchBoundaries()
   const
{
   return(d_patch_boundaries);
}


template<int DIM> inline
const IntVector<DIM>& PatchGeometry<DIM>::getRatio() const
{
   return(d_ratio_to_level_zero);
}

template<int DIM> inline
bool PatchGeometry<DIM>::intersectsPhysicalBoundary() const
{
   return(d_has_regular_boundary);
}

template<int DIM> inline
bool PatchGeometry<DIM>::getTouchesRegularBoundary() const
{
   return(d_has_regular_boundary);
}

template<int DIM> inline
bool PatchGeometry<DIM>::getTouchesPeriodicBoundary() const
{
   return(d_has_periodic_boundary);
}

template<int DIM> inline
bool
PatchGeometry<DIM>::getTouchesRegularBoundary(int axis, int upperlower) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(axis >= 0 && axis < DIM);
   TBOX_ASSERT(upperlower == 0 || upperlower == 1);
#endif

   return(d_touches_regular_bdry[axis][upperlower]);
}
                                                                                
template<int DIM> inline
bool
PatchGeometry<DIM>::getTouchesPeriodicBoundary(int axis, int upperlower) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(axis >= 0 && axis < DIM);
   TBOX_ASSERT(upperlower == 0 || upperlower == 1);
#endif

   return(d_touches_periodic_bdry[axis][upperlower]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getNodeBoundaries() const
{
   return(d_patch_boundaries[DIM-1]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getEdgeBoundaries() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(DIM >= 2);
#endif

   // The "funny" indexing prevents a warning when compiling for 
   // DIM < 2.  This code is only reached if DIM >= 2 when 
   // executing.
   return(d_patch_boundaries[DIM < 2 ? 0 : DIM-2]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getFaceBoundaries() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(DIM >= 3);
#endif

   // The "funny" indexing prevents a warning when compiling for 
   // DIM < 3.  This code is only reached if DIM >= 3 when 
   // executing.
   return(d_patch_boundaries[DIM < 3 ? 0 : DIM-3]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getCodimensionBoundaries(const int codim) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((codim > 0) && (codim <= DIM));
#endif
   return(d_patch_boundaries[codim-1]);
}

// Deprecated methods
#if (INCLUDE_DEPRECATED <= 2)
template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >* PatchGeometry<DIM>::getPatchBoundary()
   const
{
   return getPatchBoundaries();
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getNodeBoundary() const
{
   return(getNodeBoundaries());
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getEdgeBoundary() const
{
   return(getEdgeBoundaries());
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getFaceBoundary() const
{
   return(getFaceBoundaries());
}


template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getCodimensionBoundary(const int codim) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((codim > 0) && (codim <= DIM));
#endif
   return(d_patch_boundaries[codim-1]);
}

#endif // Deprecated

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
