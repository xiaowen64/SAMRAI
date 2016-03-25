//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/index/IndexDataFactory.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: hier::Patch data factory for irregularly indexed patch data
//

#ifndef included_pdat_IndexDataFactory
#define included_pdat_IndexDataFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class IndexDataFactory<DIM> is the factory class used to allocate
 * new instances of IndexData<DIM> objects.  It is a subclass of the
 * hier::PatchDataFactory<DIM> base class.
 *
 * @see pdat::IndexData
 * @see pdat::IndexVariable
 * @see hier::PatchDataFactory
 */

template<int DIM, class TYPE>
class IndexDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /**
    * The default constructor for the IndexDataFactory<DIM> class.
    * The ghost cell width argument gives the default width for all
    * irregular data objects created with this factory.
    */
   IndexDataFactory(const hier::IntVector<DIM>& ghosts);

   /**
    * Virtual destructor for the irregular data factory class.
    */
   virtual ~IndexDataFactory<DIM,TYPE>();

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
    * Virtual factory function to allocate a concrete irregular data object.
    * The default information about the object (e.g., ghost cell width) is
    * taken from the factory.  If no memory pool is provided, the allocation
    * routine assumes some default memory pool.
    */
   virtual tbox::Pointer< hier::PatchData<DIM> >
   allocate(const hier::Box<DIM>& box,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const;

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
   virtual tbox::Pointer< hier::BoxGeometry<DIM> >
   getBoxGeometry(const hier::Box<DIM>& box) const;

   /**
    * Calculate the amount of memory needed to store the irregular data
    * object, including object data but not dynamically allocated data.
    * Because the irregular data list can grow and shrink, it would be
    * impossible to estimate the necessary amount of memory.  Instead,
    * dynamic data is allocated via the standard new/free mechanisms.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /**
    * Return a boolean true value indicating that the index data quantities will always 
    * be treated as though fine values represent them on coarse-fine interfaces.
    * See the IndexVariable<DIM> class header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return false since the index data index space matches the cell-centered
    * index space for AMR patches.  Thus, index data does not live on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return false;}

   /**
    * Return whether it is valid to copy this IndexDataFactory to the 
    * supplied destination patch data factory.  It will return true if 
    * dst_pdf is an IndexDataFactory of the same type and dimension,
    * false otherwise.
    */
   bool validCopyTo(
      const tbox::Pointer< hier::PatchDataFactory<DIM> >& dst_pdf) const;   


private:

   IndexDataFactory(const IndexDataFactory<DIM,TYPE>&);
   void operator=(const IndexDataFactory<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "IndexDataFactory.C"
#endif
