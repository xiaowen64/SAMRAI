//
// File:	PatchGeometry.I
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 729 $
// Modified:	$Date: 2005-11-11 13:49:10 -0800 (Fri, 11 Nov 2005) $
// Description: Base class for geometry management on patches
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >* PatchGeometry<DIM>::getPatchBoundary()
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
   assert(axis >= 0 && axis < DIM);
   assert(upperlower == 0 || upperlower == 1);
#endif

   return(d_touches_regular_bdry[axis][upperlower]);
}
                                                                                
template<int DIM> inline
bool
PatchGeometry<DIM>::getTouchesPeriodicBoundary(int axis, int upperlower) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(axis >= 0 && axis < DIM);
   assert(upperlower == 0 || upperlower == 1);
#endif

   return(d_touches_periodic_bdry[axis][upperlower]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getNodeBoundary() const
{
   return(d_patch_boundaries[DIM-1]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getEdgeBoundary() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(DIM >= 2);
#endif

   // The "funny" indexing prevents a warning when compiling for 
   // DIM < 2.  This code is only reached if DIM >= 2 when 
   // executing.
   return(d_patch_boundaries[DIM < 2 ? 0 : DIM-2]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getFaceBoundary() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(DIM >= 3);
#endif

   // The "funny" indexing prevents a warning when compiling for 
   // DIM < 3.  This code is only reached if DIM >= 3 when 
   // executing.
   return(d_patch_boundaries[DIM < 3 ? 0 : DIM-3]);
}

template<int DIM> inline
const tbox::Array< BoundaryBox<DIM> >&
PatchGeometry<DIM>::getCodimensionBoundary(const int codim) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((codim > 0) && (codim <= DIM));
#endif
   return(d_patch_boundaries[codim-1]);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
