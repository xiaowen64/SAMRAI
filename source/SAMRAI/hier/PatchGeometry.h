/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches
 *
 ************************************************************************/

#ifndef included_hier_PatchGeometry
#define included_hier_PatchGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/LocalId.h"
#include "SAMRAI/hier/PatchBoundaries.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/List.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/**
 * Class PatchGeometry is the base class for geometry classes that
 * manage index space and mesh increment information on individual patches.
 * Patch geometry information is used for setting boundary conditions in
 * ghost cells and is used in the inter-level transfer operators for refining
 * or coarsening data between two patches associated with different index
 * spaces.  The boundary information for patches is actually computed by
 * the GridGeometry class.
 *
 * @see hier::BoundaryBox
 * @see hier::GridGeometry
 */

class PatchGeometry:public tbox::DescribedClass
{
public:
   /*!
    * @brief Array of 2*DIM booleans (with default constructor),
    * used to instantiate the sparse container map<LocalId,TwoDimBool>
    * (map<LocalId,bool[2*DIM]> does not work).
    */
   class TwoDimBool
   {
public:
      explicit TwoDimBool(
         const tbox::Dimension& dim);

      explicit TwoDimBool(
         const tbox::Dimension& dim,
         bool v);

      void
      setAll(
         bool v);

      bool&
      operator () (
         int dim,
         int side);

      const bool&
      operator () (
         int dim,
         int side) const;

      /**
       * Return the dimension of this object.
       */
      const tbox::Dimension&
      getDim() const;

      friend class::std::map<LocalId, SAMRAI::hier::PatchGeometry::TwoDimBool>;

private:
      /*
       * Needed by brain dead STL. Don't use for other purposes.
       */
      TwoDimBool();

      const tbox::Dimension d_dim;
      bool d_data[2 * tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   };

   /**
    * The default constructor for the patch geometry base class.
    */
   explicit PatchGeometry(
      const IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const TwoDimBool& touches_periodic_bdry);

   /**
    * The virtual destructor for the patch geometry base class.
    */
   virtual ~PatchGeometry();

   /**
    * Return const reference to patch boundary information.
    */
   const tbox::Array<tbox::Array<BoundaryBox> >
   getPatchBoundaries() const;

   /*!
    * @brief Set the boundary box arrays for this patch geometry.
    *
    * An array of length DIM of tbox::Array< BoundaryBox > is passed
    * in to be stored as the boundary boxes for this patch geometry.
    *
    * @param bdry The array of BoundaryBox arrays.
    */
   void
   setBoundaryBoxesOnPatch(
      const tbox::Array<tbox::Array<BoundaryBox> > bdry);

   /**
    * Return const reference to ratio to level zero index space.
    */
   const IntVector&
   getRatio() const;

   /**
    * Return a boolean value indicating whether the patch boundary
    * intersects the physical domain boundary in a non-periodic
    * direction.  In other words, the return value is true when the
    * patch has non-empty boundary boxes that lie outside the physical
    * domain.  Otherwise, the return value is false.  Note that when
    * a patch touches the "boundary" of the physical domain in a periodic
    * direction, there are no boundary boxes to fill; the data is filled
    * from the proper region of the domain interior in the periodic direction.
    */
   bool
   intersectsPhysicalBoundary() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch at a single point (i.e., 0-dim intersection
    * between cells in patch and cells in boundary box).
    */
   const tbox::Array<BoundaryBox>&
   getNodeBoundaries() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch along a 1-dim edge (i.e., 1-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * When assertion checking is active, this routine throws an assertion
    * when DIM < 2.
    */
   const tbox::Array<BoundaryBox>&
   getEdgeBoundaries() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch along a 2-dim face (i.e., 2-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * When assertion checking is active, this routine throws an assertion
    * when DIM < 3.
    */
   const tbox::Array<BoundaryBox>&
   getFaceBoundaries() const;

   /**
    * Return array of boundary box components for patch each of which
    * intersects the patch as a (DIM - codim)-dimensional object.
    * That is,
    *
    * if DIM == 1: (co(dim == tbox::Dimension(1))) => same components as getNodeBoundaries.
    *
    * if DIM == 2, (co(dim == tbox::Dimension(1))) => same components as getEdgeBoundaries.
    *              (co(dim == tbox::Dimension(2))) => same components as getNodeBoundaries.
    *
    * if DIM == 3, (co(dim == tbox::Dimension(1))) => same components as getFaceBoundaries.
    *              (co(dim == tbox::Dimension(2))) => same components as getEdgeBoundaries.
    *              (co(dim == tbox::Dimension(3))) => same components as getNodeBoundaries.
    *
    * When assertion checking is active, this routine throws an assertion
    * when codim < 0 or codim > DIM.
    */
   const tbox::Array<BoundaryBox>&
   getCodimensionBoundaries(
      const int codim) const;

   /**
    * Set the array of boundary box components of the given codimension
    * for a patch.
    */
   void
   setCodimensionBoundaries(
      const tbox::Array<BoundaryBox>& bdry_boxes,
      const int codim);

   /*!
    * @brief Compute a box outside a physical domain that needs to be filled.
    *
    * The patch box will be grown by the given ghost cell width and
    * then intersected with the boundary box.  The resulting intersection
    * will be grown to the needed ghost cell width in the direction
    * normal to the boundary.
    *
    * @param bbox BoundaryBox representing location and type of boundary
    * @param patch_box The box for the patch where data is being filled
    * @param gcw ghost cell width to fill
    */
   Box
   getBoundaryFillBox(
      const BoundaryBox& bbox,
      const Box& patch_box,
      const IntVector& gcw) const;

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any non-periodic physical boundary
    */
   bool
   getTouchesRegularBoundary() const;

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any periodic boundary
    */
   bool
   getTouchesPeriodicBoundary() const;

   /*!
    * @brief Query whether patch touches a specific regular boundary
    *
    * Returns true if the Patch touches a non-periodic physical boundary
    * on the side of the Patch specified in the argument list.  The side
    * is specified by an axis direction and a flag specified the upper or
    * lower side.
    *
    * @param axis       Axis direction normal to the side being checked
    * @param upperlower Flag should be 0 if checking the lower side in the
    *                   axis direction, or 1 if checking the upper side.
    */
   bool
   getTouchesRegularBoundary(
      int axis,
      int upperlower) const;

   /**
    * Print object data to the specified output stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

private:
   const tbox::Dimension d_dim;

   bool d_has_regular_boundary;
   bool d_has_periodic_boundary;
   IntVector d_ratio_to_level_zero;
   PatchBoundaries d_patch_boundaries;

   TwoDimBool d_touches_regular_bdry;

   tbox::List<IntVector> d_periodic_shifts;
};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/PatchGeometry.I"
#endif
#endif
