//
// File:	CellDataFactory.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Factory class for creating cell data objects
//

#ifndef included_pdat_CellDataFactory
#define included_pdat_CellDataFactory

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
 * Class CellDataFactory is a factory class used to allocate new
 * instances of CellData objects.  It is a subclass of the patch
 * data factory class and cell data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::CellData
 * @see pdat::PatchDataFactory
 */

template<int DIM, class TYPE>
class CellDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The default constructor for the cell data factory class.  The ghost
    * cell width and depth (number of components) arguments give the defaults
    * for all cell data objects created with this factory.
    */
   CellDataFactory(int depth, 
                         const hier::IntVector<DIM>& ghosts);

   /**
    * Virtual destructor for the cell data factory class.
    */
   virtual ~CellDataFactory<DIM,TYPE>();

   /**
    * Virtual function to clone the patch data factory .  This will return
    * a new instantiation of the factory with the same properties (e.g., same
    * type and ghost cell width).  The properties of the cloned factory can
    * then be changed without modifying the original.
    */
   virtual tbox::Pointer< hier::PatchDataFactory<DIM> > cloneFactory();

   /**
    * Virtual factory function to allocate a concrete cell data object.
    * The default information about the object (e.g., ghost cell width)
    * is taken from the factory.  If no memory pool is provided, then
    * the allocation routine assumes some default memory pool.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> > allocate(
      const hier::Box<DIM>& box,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    */

   virtual tbox::Pointer< hier::BoxGeometry<DIM> > getBoxGeometry(
      const hier::Box<DIM>& box) const;

   /**
    * Get the default ghost cell width.  This is the ghost cell width that
    * will be used in the instantiation of concrete cell data objects.
    */
   virtual const hier::IntVector<DIM>& getDefaultGhostCellWidth() const;

   /**
    * Set the default ghost cell width.  This is the ghost cell width that
    * will be used in the instantiation of concrete cell data instances.
    */
   virtual void setDefaultGhostCellWidth(const hier::IntVector<DIM>& ghosts);

   /**
    * Get the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of cell data objects.
    */
   int getDefaultDepth() const;

   /**
    * Set the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of cell data objects.
    */
   void setDefaultDepth(const int depth);

   /**
    * Calculate the amount of memory needed to store the cell data object,
    * including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean true value indicating that the cell data quantities will always
    * be treated as though fine values represent them on coarse-fine interfaces.
    * See the CellVariable<DIM> class header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return false since the cell data index space matches the cell-centered
    * index space for AMR patches.  Thus, cell data does not live on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return false;}

   /**
    * Return whether it is valid to copy this CellDataFactory to the 
    * supplied destination patch data factory. It will return true if
    * dst_pdf is a CellDataFactory, false otherwise.
    */
   bool validCopyTo(
      const tbox::Pointer< hier::PatchDataFactory<DIM> >& dst_pdf) const;   


private:
   int d_depth;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellDataFactory.C"
#endif
