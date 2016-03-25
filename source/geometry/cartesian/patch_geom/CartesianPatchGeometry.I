//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/geometry/cartesian/patch_geom/CartesianPatchGeometry.I $
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Simple Cartesian grid geometry for an AMR hierarchy.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace geom {

template<int DIM> inline 
const double* CartesianPatchGeometry<DIM>::getDx() const
{
   return(d_dx);
} 

template<int DIM> inline 
const double* CartesianPatchGeometry<DIM>::getXLower() const
{
   return(d_x_lo);
}
 
template<int DIM> inline 
const double* CartesianPatchGeometry<DIM>::getXUpper() const
{
   return(d_x_up);
}
 
#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
