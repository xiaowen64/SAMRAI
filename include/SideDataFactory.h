//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/side/SideDataFactory.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description: Factory class for creating side data objects
//

#ifndef included_pdat_SideDataFactory
#define included_pdat_SideDataFactory

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxGeometry
#include "BoxGeometry.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_PatchDataFactory
#include "PatchDataFactory.h"
#endif
#ifndef included_pdat_MultiblockSideDataTranslator
#include "MultiblockSideDataTranslator.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class SideDataFactory is a factory class used to allocate new
 * instances of SideData objects.  It is a subclass of the patch
 * data factory class and side data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * Note that it is possible to create a side data factory to allocate
 * and manage data for cell sides associated with a single coordinate 
 * direction only.  See the constructor for more information.
 *
 * @see pdat::SideData
 * @see pdat::PatchDataFactory
 */

template<int DIM, class TYPE>
class SideDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The default constructor for the side data factory class.  The ghost cell 
    * width, depth (number of components), and fine boundary representation arguments 
    * give the defaults for all edge data objects created with this factory.
    * Also, the default data allocation scheme is to generate storage for sides 
    * in all coordinate directions (default integer vector of all 1's).  To 
    * use this factory to manage side data objects for sides associated
    * with a single direction only, provide the directions vector argument.
    * A zero entry indicates that data for that direction is not wanted.
    * Otherwise, data will be created for that direction.  See the 
    * SideVariable<DIM> class header file for more information.
    */
   SideDataFactory(int depth, 
                         const hier::IntVector<DIM>& ghosts,
                         bool fine_boundary_represents_var,
                         const hier::IntVector<DIM>& directions = 
                                                hier::IntVector<DIM>(1));

   /**
    * Virtual destructor for the side data factory class.
    */
   virtual ~SideDataFactory<DIM,TYPE>();

   /**
    * @brief Abstract virtual function to clone a patch data factory.

    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts default ghost cell width for concrete classes created from
    * the factory.
    */
   virtual tbox::Pointer< hier::PatchDataFactory<DIM> > cloneFactory(const hier::IntVector<DIM>& ghosts);

   /**
    * Virtual factory function to allocate a concrete side data object.
    * The default information about the object (e.g., ghost cell width)
    * is taken from the factory.  If no memory pool is provided, then
    * the allocation routine assumes some default memory pool.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> > allocate(
      const hier::Box<DIM>& box,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /**
    * Virtual factory function to allocate a concrete cell data object.
    * Same as above function, except passes in a patch instead of a box.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> > allocate(
      const hier::Patch<DIM>& patch,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    */
   virtual tbox::Pointer< hier::BoxGeometry<DIM> > getBoxGeometry(
      const hier::Box<DIM>& box) const;

   /**
    * Get the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of side data objects.
    */
   int getDefaultDepth() const;

   /**
    * Set the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of side data objects.
    */
   void setDefaultDepth(const int depth);

   /**
    * Return constant reference to vector describing which coordinate
    * directions have data associated with this side data object.
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction.  A non-zero
    * value indicates that a valid data array is maintained for that
    * coordinate direction.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

   /**
    * Calculate the amount of memory needed to store the side data object,
    * including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean value indicating how data for the side quantity will be treated
    * on coarse-fine interfaces.  This value is passed into the constructor.  See
    * the FaceVariable<DIM> class header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return d_fine_boundary_represents_var;}

   /**
    * Return true since the side data index space extends beyond the interior of
    * patches.  That is, side data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

   /**
    * Return whether it is valid to copy this SideDataFactory to the 
    * supplied destination patch data factory.  It will return true if 
    * dst_pdf is SideDataFactory or OutersideDataFactory, false otherwise.
    */
   bool validCopyTo(
      const tbox::Pointer< hier::PatchDataFactory<DIM> >& dst_pdf) const;   

   /**
    * Return pointer to a multiblock data translator
    */
   hier::MultiblockDataTranslator<DIM>* getMultiblockDataTranslator();

private:
   int d_depth;
   bool d_fine_boundary_represents_var;
   hier::IntVector<DIM> d_directions;

   MultiblockSideDataTranslator<DIM,TYPE>* d_mb_trans;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "SideDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideDataFactory.C"
#endif
