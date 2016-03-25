//
// File:	PatchGeometry.C
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Base class for geometry management on patches
//

#ifndef included_hier_PatchGeometry_C
#define included_hier_PatchGeometry_C

#include "PatchGeometry.h"

#include "BoundaryLookupTable.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "PatchGeometry.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchGeometry<DIM>::PatchGeometry(
   const IntVector<DIM>& ratio_to_level_zero, 
   const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
   const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must me nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < DIM; i++) {
      assert( ratio_to_level_zero(i) != 0 );
   }
   if (DIM > 1) {
      for (i = 0; i < DIM; i++) {
	 assert( (ratio_to_level_zero(i)*ratio_to_level_zero((i+1)%DIM) > 0)
		 || (ratio_to_level_zero(i) == 1)
		 || (ratio_to_level_zero((i+1)%DIM) == 1) );
      }
   }

   assert(touches_regular_bdry.size() == DIM);
   assert(touches_periodic_bdry.size() == DIM);

   for (int i = 0; i < DIM; i++) {
      assert(touches_regular_bdry[i].size() == 2); 
      assert(touches_periodic_bdry[i].size() == 2); 
   }
#endif

   d_has_regular_boundary = false;
   d_has_periodic_boundary = false;
   d_ratio_to_level_zero = ratio_to_level_zero;

   for (int axis = 0; axis < DIM; axis++) {
      for (int dir = 0; dir < 2; dir++) {
         d_touches_regular_bdry[axis][dir] = touches_regular_bdry[axis][dir];
         d_touches_periodic_bdry[axis][dir] = touches_periodic_bdry[axis][dir];

         if (d_touches_regular_bdry[axis][dir]) {
            d_has_regular_boundary = true;
         }
         if (d_touches_periodic_bdry[axis][dir]) {
            d_has_periodic_boundary = true;
         }
      }
   }
}

template<int DIM>  PatchGeometry<DIM>::~PatchGeometry()
{
}

template<int DIM> Box<DIM>
PatchGeometry<DIM>::getBoundaryFillBox(const BoundaryBox<DIM>& bbox,
                                        const Box<DIM>& patch_box,
                                        const IntVector<DIM>& gcw) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < DIM; i++) {
      assert( gcw(i) >= 0 );
   }
#endif

   Box<DIM> tmp_box(patch_box);
   tmp_box.grow(gcw);
   Box<DIM> fill_box(bbox.getBox() * tmp_box);

   int bdry_type = bbox.getBoundaryType();
   int location_index = bbox.getLocationIndex();

   // Get the singleton class lookup table 
   const BoundaryLookupTable<DIM>* blut;
   blut = BoundaryLookupTable<DIM>::getLookupTable();

#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Array<int> &location_index_max = blut->getMaxLocationIndices();
   assert(bdry_type > 0);
   assert(bdry_type <= DIM);
   assert(location_index >= 0);
#endif

   if (!fill_box.empty()) {

      // Loop over codimension (a.k.a. boundary type)
      for (int codim = 1; codim <= DIM; codim++) {

	 // When we get a match on the boundary type
         if (bdry_type == codim) {

#ifdef DEBUG_CHECK_ASSERTIONS
            assert(location_index < location_index_max[codim-1]);
#endif
            // Get the directions involved in this boundary type from the 
	    // lookup table.
	    const tbox::Array<int> &dir =
               blut->getDirections(location_index, codim);

	    // For each direction, identify this as an upper or lower boundary.
            for (int i = 0; i < codim; i++) {
	       if (blut->isUpper(location_index, codim, i)) {
                  fill_box.growUpper(dir[i], gcw(dir[i])-1);
	       } else {
                  fill_box.growLower(dir[i], gcw(dir[i])-1);
	       }
            }

	    // We've found boundary type, so break out of the loop.
            break;
         }
      }
   }

   return (fill_box);
}


template<int DIM> void
PatchGeometry<DIM>::setCodimensionBoundary(
   const tbox::Array< BoundaryBox<DIM> >& bdry_boxes,
   int codim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < bdry_boxes.size(); i++) {
      assert(bdry_boxes[i].getBoundaryType() == codim);
   }
   assert(codim <= DIM);
   assert(codim > 0);
#endif
                                                                                
   d_patch_boundaries[codim-1].resizeArray(bdry_boxes.size());
                                                                                
   for (int b = 0; b < bdry_boxes.size(); b++) {
      d_patch_boundaries[codim-1][b] = bdry_boxes[b];
   }
}

template<int DIM> void
PatchGeometry<DIM>::setBoundaryBoxesOnPatch(
   const tbox::Array< BoundaryBox<DIM> > bdry[DIM])
{
   for (int i = 0; i < DIM; i++) {
      setCodimensionBoundary(bdry[i], i+1);
   }
}

template<int DIM> void PatchGeometry<DIM>::printClassData(ostream& stream) const
{
   stream << "\nPatchGeometry<DIM>::printClassData..." << endl;
   stream << "Ratio to level zero = " << d_ratio_to_level_zero << endl;
   stream << "d_has_regular_boundary = " 
          << d_has_regular_boundary << endl;
   stream << "Boundary boxes for patch..." << endl;
   for (int d = 0; d < DIM; d++) {
      const int n = d_patch_boundaries[d].getSize();
      stream << "Boundary box array " << d << " has " << n << " boxes" << endl;
      for (int i = 0; i < n; i++) {
         stream << "box " << i << " = "
                << d_patch_boundaries[d][i].getBox() << endl;
      }
   }
}

}
}
#endif
