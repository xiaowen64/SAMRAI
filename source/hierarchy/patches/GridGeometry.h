//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/GridGeometry.h $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Base class for geometry management in AMR hierarchy
//

#ifndef included_hier_GridGeometry
#define included_hier_GridGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <std::string>
#define included_String
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM> class PatchLevel;

/**
 * Class GridGeometry<DIM> serves as the base class for SAMRAI geometry
 * classes that manage particular grid types (e.g., Cartesian, cylindrical,
 * etc.).  The grid geometry class is responsible for maintaining information
 * about the index space describing the physical domain and computing this
 * information for patches in an AMR hierarchy.  Operations performed by
 * this class include determining which patches are adjacent to the physical
 * domain boundary and computing boundary boxes for patches which decribe
 * how the patch touches the domain boundary (useful for filling ghost cell
 * data for physical boundary conditions).  Member functions that
 * manage the description of the spatial coordinated on the mesh are pure
 * virtual here and must be implemented in an appropriate subclass.
 *
 * Note that the derivation from tbox::DescribedClass is virtual.  The
 * reason for this is to avoid dynamic casting problems for smart pointers.
 * Typically, SAMRAI geometry objects inherit from tbox::Serializable
 * as well as this base geometry class.  Thus, there is usually more than
 * one class hierarchy for geometry objects.  Pointers to base objects
 * may need to be dynamically cast to derived objects in either hierarchy.
 *
 * @see hier::BoundaryBox
 */

template<int DIM> class GridGeometry : public virtual tbox::DescribedClass
{
public:
   /**
    * Constructor for GridGeometry.
    */
   GridGeometry(const std::string &object_name);

   /**
    * Destructor for GridGeometry.
    */
   virtual ~GridGeometry();

   //! @{
   /*!
    * @name Functions for computing boundary boxes
    */
   /*!
    * @brief Determine for every patch on a level if it touches a regular
    * physical boundary or a periodic boundary.
    *
    * This routine loops through all of the patches on the given level
    * and determines which kinds of boundaries each patch touches.  The
    * 3-dimensional boolean arrays are set to store for each path whether
    * it touches a regular boundary, a periodic boundary, both, or
    * neither.
    *
    * The array arguments should be uninitialized when they are passed into
    * this function.
    *
    * @param touches_regular_bdry Array to store which patches touch
    *                             non-periodic boundaries.
    * @param touches_periodic_bdry Array to store which patches touch
    *                              periodic boundaries.
    * @param level containing the patches to be checked
    * @param periodic_shift periodic shift for the level (see getPeriodicShift)
    * @param domain Physical domain (at the same level of refinement as level)
    */
   void findPatchesTouchingBoundaries(
      tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_regular_bdry,
      tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_periodic_bdry,
      const PatchLevel<DIM>& level,
      const IntVector<DIM>& periodic_shift,
      const BoxArray<DIM>& domain) const;

   /*!
    * @brief Pass the arrays holding the boundary information to be stored 
    * in the concrete geometry classes, and construct boundary boxes if
    * required
    *
    * This routine will pass the arrays containing the information about
    * which patches touch which boundaries to the concrete grid geometry
    * class.  Also, if defer_boundary_box_creation is false, this routine
    * will call a routine to construct all of the boundary boxes for the
    * patches on level and will set them in the patch geometry for each patch.
    *
    * @param level containing the patches to be checked.
    * @param ratio_to_level_zero ratio to the coarsest level.
    * @param touches_regular_bdry Array storing which patches touch
    *                             non-periodic boundaries.
    * @param touches_periodic_bdry Array storing which patches touch
    *                              periodic boundaries.
    * @param defer_boundary_box_creation Boundary boxes will be created here
    *                                    if false, not if true.
    */
   virtual void setGeometryOnPatches(
      hier::PatchLevel<DIM>& level,
      const hier::IntVector<DIM>& ratio_to_level_zero,
      tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_regular_bdry,
      tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_periodic_bdry,
      bool defer_boundary_box_creation);

   /*!
    * @brief Construct the boundary boxes for each patch and set them on
    * the patch geometries.
    *
    * This routine constructs the boundary boxes for every patch in the level.
    * Once constructed, the boundary boxes are set on each patch's
    * PatchGeometry object.
    *
    * @param level The level for which boundary boxes are constructed.
    */
   void setBoundaryBoxes(hier::PatchLevel<DIM>& level);

   //@}

   /**
    * Compute the valid periodic shifts for each patch on a level.  The
    * shifts array will store a list of IntVectors for each patch.
    * Each list will contain the valid possible periodic shifts for each
    * particular patch.  If there are no periodic boundary conditions or a
    * patch does not touch a periodic boundary, the list for a patch will be
    * empty.
    *
    * The patch geometry object for each patch on the level must be properly
    * initialized before calling this routine.  This is typically done in 
    * the patch level constructor.
    *
    * When assertion checking is active, the array of shifts in the
    * argument list must have the same length as the number of patches
    * on the level.
    */
   void computeShiftsForLevel(
      tbox::Array< tbox::List< IntVector<DIM> > >& shifts,
      const PatchLevel<DIM>& level,
      const BoxArray<DIM>& physical_domain) const;

   /**
    * Compute physical domain box array describing the index space
    * of the physical domain managed by this geometry object.  If any entry of
    * ratio vector is negative, the index space is coarsened with respect to
    * the physical domain description.  Otherwise, the index space is refined.
    */
   void computePhysicalDomain(
      BoxArray<DIM>& domain,
      const IntVector<DIM>& ratio_to_level_zero) const;

   /**
    * Set physical domain to input box array and determine whether
    * domain is a single box.
    */
   void setPhysicalDomain(const BoxArray<DIM>& domain);

   /**
    * Return const reference to physical domain description for level 0.
    */
   const BoxArray<DIM>& getPhysicalDomain() const;

   /**
    * Return boolean value indicating whether the physical domain can be
    * represented as a single box.
    */
   bool getDomainIsSingleBox() const;

   /**
    * Print object data to the specified output stream.
    */
   virtual void printClassData(std::ostream& stream) const;

   /**
    * Initialize the periodic shift on the coarsest level.  The IntVector
    * argument should be set to 1 for periodic directions and 0 for all
    * other directions.  The shift will be calculated to be the number of
    * cells in the periodic direction and zero in all other directions.
    */
   void initializePeriodicShift(const IntVector<DIM>& directions);

   /**
    * Return IntVector<DIM> containing the periodic shift in each direction
    * for a domain represented by a refinement of the reference physical
    * domain (i.e. level zero) by the given ratio vector.  tbox::Array entries
    * will be zero for non-periodic directions.  By default (i.e., when
    * no argument is passed, the function returns the periodic shift for
    * level zero in the hierarchy.
    */
   IntVector<DIM> getPeriodicShift(
      const IntVector<DIM>& ratio_to_level_zero = IntVector<DIM>(1)) const;

   /*!
    * @brief Compute the maximum ghost width of all of the
    * components associated with the patch descriptor.
    *
    * Calculates the maximum ghost width for all the variables associated
    * with the patch descriptor.  This must only be called after all of the
    * variables have been registered with the VariableDatabase.  If a
    * variable is added that changes the maximum ghost width, then an
    * assertion failure will result.
    */
   IntVector<DIM> computeMaxGhostWidth(
      tbox::Pointer< PatchDescriptor<DIM> > descriptor);

   /**
    * Pure virtual function to create and return a pointer to a refined
    * version of this grid geometry object.
    */
   virtual tbox::Pointer<hier::GridGeometry<DIM> > makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector<DIM>& refine_ratio,
      bool register_for_restart) const = 0;

   /**
    * Pure virtual function to create and return a pointer to a coarsened
    * version of this grid geometry object.
    */
   virtual tbox::Pointer<hier::GridGeometry<DIM> > makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector<DIM>& coarsen_ratio,
      bool register_for_restart) const = 0;

   /*
    * Pure virtual function to compute grid data for patch and assign new
    * concrete patch geometry object to patch.
    */
   virtual void setGeometryDataOnPatch(
      hier::Patch<DIM>& patch,
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry) const = 0;

   /*!
    * @brief Compute boundary boxes for each patch in patch level and
    * assign them to the array of boundary box arrays, assumed to be
    * of length DIM * (num patches).
    *
    * The DIM arrays of boundary boxes for each patch will be stored in groups
    * of DIM. For example, in 3d with n patches on the level, the array
    * For example, in 3d with n patches on the level, the array
    * of boundary box arrays will be ordered as follows:
    *
    * @verbatim
    * (patch 0 face array, patch 0 edge array, patch 0 node array,
    *  patch 1 face array, patch 1 edge array, patch 1 node array, . . . ,
    *  patch n-1 face array, patch n-1 edge array, patch n-1 node array)
    * @endverbatim
    *
    * The optional argument do_all_patches defaults to false, in which case
    * the boundary box computation is executed only on patches that touch a
    * non-periodic boundary.  When this routine is called during patch
    * level construction to describe a physical boundary, it is known that
    * only patches that touch a non-periodic boundary will have non-empty
    * sets of boundary boxes, so for efficiency's sake the boundary box
    * box computation is supressed for all other patches.  When this
    * routine is called to create boundary boxes that describe a
    * coarse-fine boundary, the computation must occur for every patch, so
    * do_all_patches mush be set to true.
    *
    * @param boundaries output boundary description
    * @param level level on which to generate boundaries
    * @param periodic_shift periodic shift for the level (see getPeriodicShift)
    * @param ghost_width ghost width to compute geometry for
    * @param domain Physical domain (in index space of level) for computing
    *               boundary boxes.
    * @param do_all_patches Execute boundary box computation on all patches,
    *                       even those known to not touch a boundary
    */
   void computeBoundaryBoxesOnLevel(
      tbox::Array< BoundaryBox<DIM> > boundaries[],
      const PatchLevel<DIM>& level,
      const IntVector<DIM>& periodic_shift,
      const IntVector<DIM>& ghost_width,
      const BoxArray<DIM>& domain,
      bool do_all_patches = false) const;

   /*!
    * @brief Compute boundary boxes for patch
    *
    * Decompose patch boundary region into pieces depending on spatial dim.
    * Boxes are extended along the boundary to the edge of the ghost layer
    * if necessary.
    */
   void getBoundaryBoxes(
      tbox::Array< BoundaryBox<DIM> > boundaries[DIM],
      const Box<DIM>& box,
      const BoxArray<DIM>& domain_boxes,
      const IntVector<DIM>& ghosts,
      const IntVector<DIM> &periodic_shift) const;

private:

   /*!
    * @brief Check that the domain is valid for periodic boundary conditions
    */
   bool checkPeriodicValidity(const BoxArray<DIM>& domain);

   /*!
    * @brief Check on each BoundaryBox when it is created.
    *
    * This is a check performed on each BoundaryBox when it is created.
    * It returns true when a BoundaryBox has a width of 1 in at least
    * one direction, is adjacent to the patch boundary (possible extended
    * into the patch's ghost region) and is outside the physical domain.
    */
   bool checkBoundaryBox(
      const BoundaryBox<DIM>& boundary_box,
      const Patch<DIM>& patch,
      const BoxArray<DIM>& domain,
      const int num_per_dirs,
      const IntVector<DIM> &max_data_ghost_width ) const;

   /*!
    * @brief Find if a box is on a periodic boundary and compute its shifts.
    *
    * If box is located on a periodic boundary, all of its possible shifts
    * will be computed and stored in shifts.  If box is not on a periodic
    * boundary, shifts will be an empty list.
    */
   void computeShiftsForPatch(
      tbox::List< IntVector<DIM> >& shifts,
      const Box<DIM>& box,
      const BoxArray<DIM>& domain,
      const IntVector<DIM>& periodic_shift) const;

   /*!
    * Object name used for error reporting purposes.
    */
   std::string d_object_name;

   /*!
    * Box array defining computational domain on coarsest level
    * and boolean flag that is true when domain is a single box.
    */
   BoxArray<DIM> d_physical_domain;
   bool d_domain_is_single_box;

   /*!
    * Integer array vector describing periodic shift coarsest level.
    * An entry of zero means direction is not periodic.
    */
   IntVector<DIM> d_periodic_shift;

   /*!
    * Current maximum ghost cell width over all patch data objects
    * known to the patch descriptor.  This is used to compute
    * boundary boxes.
    */
   IntVector<DIM> d_max_data_ghost_width;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "GridGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "GridGeometry.C"
#endif
