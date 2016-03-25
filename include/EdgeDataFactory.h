//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/edge/EdgeDataFactory.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Factory class for creating edge data objects
//

#ifndef included_pdat_EdgeDataFactory
#define included_pdat_EdgeDataFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "IntVector.h"
#include "PatchDataFactory.h"
#include "MultiblockEdgeDataTranslator.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class EdgeDataFactory is a factory class used to allocate new
 * instances of EdgeData objects.  It is a subclass of the patch
 * data factory class and edge data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::EdgeData
 * @see pdat::PatchDataFactory
 */

template<int DIM, class TYPE>
class EdgeDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The constructor for the edge data factory class.  The ghost cell width, depth 
    * (number of components), and fine boundary representation arguments give the 
    * defaults for all edge data objects created with this factory.  See
    * the EdgeVariable<DIM> class header file for more information.
    */
   EdgeDataFactory(int depth, 
                   const hier::IntVector<DIM>& ghosts,
                   bool fine_boundary_represents_var);

   /**
    * Virtual destructor for the edge data factory class.
    */
   virtual ~EdgeDataFactory<DIM,TYPE>();

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
    * Virtual factory function to allocate a concrete edge data object.
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
    * depth that will be used in the instantiation of edge data objects.
    */
   int getDefaultDepth() const;

   /**
    * Set the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of edge data objects.
    */
   void setDefaultDepth(const int depth);

   /**
    * Calculate the amount of memory needed to store the edge data object,
    * including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean value indicating how data for the edge quantity will be treated
    * on coarse-fine interfaces.  This value is passed into the constructor.  See 
    * the EdgeVariable<DIM> class header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return d_fine_boundary_represents_var;}

   /**
    * Return true since the edge data index space extends beyond the interior of
    * patches.  That is, edge data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

   /**
    * Return whether it is valid to copy this EdgeDataFactory to the 
    * supplied destination patch data factory.  It will return true if 
    * dst_pdf is EdgeDataFactory or OuteredgeDataFactory, false otherwise.
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

   MultiblockEdgeDataTranslator<DIM,TYPE>* d_mb_trans;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "EdgeDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeDataFactory.C"
#endif
