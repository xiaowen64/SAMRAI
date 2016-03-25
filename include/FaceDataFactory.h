//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/face/FaceDataFactory.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Factory class for creating face data objects
//

#ifndef included_pdat_FaceDataFactory
#define included_pdat_FaceDataFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "IntVector.h"
#include "PatchDataFactory.h"
#include "MultiblockFaceDataTranslator.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class FaceDataFactory is a factory class used to allocate new
 * instances of FaceData objects.  It is a subclass of the patch
 * data factory class and face data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::FaceData
 * @see pdat::PatchDataFactory
 */

template<int DIM, class TYPE>
class FaceDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The constructor for the face data factory class.  The ghost cell width, depth
    * (number of components), and fine boundary representation arguments give the
    * defaults for all edge data objects created with this factory.  See
    * the FaceVariable<DIM> class header file for more information.
    */
   FaceDataFactory(int depth, 
                         const hier::IntVector<DIM>& ghosts,
                         bool fine_boundary_represents_var);

   /**
    * Virtual destructor for the face data factory class.
    */
   virtual ~FaceDataFactory<DIM,TYPE>();

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
    * Virtual factory function to allocate a concrete face data object.
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
    * depth that will be used in the instantiation of face data objects.
    */
   int getDefaultDepth() const;

   /**
    * Set the default depth (number of components).  This is the default
    * depth that will be used in the instantiation of face data objects.
    */
   void setDefaultDepth(const int depth);

   /**
    * Calculate the amount of memory needed to store the face data object,
    * including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean value indicating how data for the face quantity will be treated
    * on coarse-fine interfaces.  This value is passed into the constructor.  See
    * the FaceVariable<DIM> class header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return d_fine_boundary_represents_var;}

   /**
    * Return true since the face data index space extends beyond the interior of
    * patches.  That is, face data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

   /**
    * Return whether it is valid to copy this FaceDataFactory to the 
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is FaceDataFactory and OuterfaceDataFactory, false otherwise.
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

   MultiblockFaceDataTranslator<DIM,TYPE>* d_mb_trans;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "FaceDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceDataFactory.C"
#endif
