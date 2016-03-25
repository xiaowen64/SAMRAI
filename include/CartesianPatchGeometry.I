//
// File:	CartesianPatchGeometry.I
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
