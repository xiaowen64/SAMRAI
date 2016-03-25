//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/geometry/skeleton/patch_geom/BlockPatchGeometry.C $
// Package:     SAMRAI geometry package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 878 $
// Modified:    $LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description: Grid geometry for multiblock.
//

#ifndef included_geom_BlockPatchGeometry_C
#define included_geom_BlockPatchGeometry_C

#include "BlockPatchGeometry.h"

namespace SAMRAI {
    namespace geom {

/*
*************************************************************************
*                                                                       *
* Constructor for BlockPatchGeometry.                           *
* variable.                                                             *
*                                                                       *
*************************************************************************
*/
template<int DIM>  BlockPatchGeometry<DIM>::BlockPatchGeometry(
   const hier::IntVector<DIM>& ratio_to_level_zero,
   const int block_number,
   const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
   const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry)
: hier::PatchGeometry<DIM>(ratio_to_level_zero,
                           touches_regular_bdry,
                           touches_periodic_bdry)
{
   d_block_number = block_number;
}


/*
*************************************************************************
*                                                                       *
* Destructor for BlockPatchGeometry.                            *
*                                                                       *
*************************************************************************
*/
template<int DIM>  BlockPatchGeometry<DIM>::~BlockPatchGeometry()
{
}


/*
*************************************************************************
*                                                                       *
* Print BlockPatchGeometry class data.                          *
*                                                                       *
*************************************************************************
*/
template<int DIM> void BlockPatchGeometry<DIM>::printClassData(std::ostream& os) const
{
   os << "Printing BlockPatchGeometry data: this = "
      << (BlockPatchGeometry<DIM>*)this << std::endl;
 
   hier::PatchGeometry<DIM>::printClassData(os);
}

}
}
#endif
