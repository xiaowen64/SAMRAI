//
// File:	OuteredgeDataFactory.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Factory class for creating outeredge data objects
//

#ifndef included_pdat_OuteredgeDataFactory
#define included_pdat_OuteredgeDataFactory

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

/*!
 * @brief
 * Class OuteredgeDataFactory<DIM> is a factory class used to allocate new
 * instances of OuteredgeData<DIM> objects.  It is a subclass of the patch
 * data factory class and outeredge data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::OuteredgeData
 * @see pdat::PatchDataFactory
 */

template <int DIM, class TYPE>
class OuteredgeDataFactory : public hier::PatchDataFactory<DIM>
{
public:
   /*!
    * @brief
    * The default constructor for the outeredge data factory class.
    *
    * The depth (number of components) gives the default for all of
    * the outeredge data objects created with this factory.
    */
   OuteredgeDataFactory(int depth);

   /*!
    * @brief
    * Virtual destructor for the outeredge data factory class.
    */
   virtual ~OuteredgeDataFactory<DIM,TYPE>();

   /*!
    * @brief
    * Virtual function to clone the patch data factory .
    *
    * This will return
    * a new instantiation of the factory with the same properties (e.g., same
    * type and depth).  The properties of the cloned factory can then be
    * changed without modifying the original.
    */
   virtual tbox::Pointer<hier::PatchDataFactory<DIM> > cloneFactory();

   /*!
    * @brief
    * Virtual factory function to allocate a concrete outeredge data object.
    *
    * The default information about the object (e.g., depth) is taken from
    * the factory.  If no memory pool is provided, then the allocation routine
    * assumes some default memory pool.
    */
   virtual tbox::Pointer<hier::PatchData<DIM> > allocate(
      const hier::Box<DIM>& box,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL)) const;

   /*!
    * @brief
    * Allocate the box geometry object associated with the patch data.
    *
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    */
   virtual tbox::Pointer<hier::BoxGeometry<DIM> >
   getBoxGeometry(const hier::Box<DIM>& box) const;

   /*!
    * @brief
    * Get the default ghost cell width.
    *
    * This is the ghost cell width that will be used in the instantiation of 
    * concrete outeredge data objects.  This is set to zero currently for 
    * outeredge data objects.
    */
   virtual const hier::IntVector<DIM>& getDefaultGhostCellWidth() const;

   /*!
    * @brief
    * Set the default ghost cell width.
    *
    * This is the ghost cell width that will be used in the instantiation of 
    * concrete outeredge data instances.  This routine should not be called 
    * with a non-zero ghost cell width and will cause an abort if the code is 
    * compiled with assertions enabled.
    */
   virtual void setDefaultGhostCellWidth(const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief
    * Get the default depth (number of components).
    *
    * This is the default depth that will be used in the instantiation of 
    * outeredge data objects.
    */
   int getDefaultDepth() const;

   /*!
    * @brief
    * Set the default depth (number of components).
    *
    * This is the default depth that will be used in the instantiation 
    * of outeredge data objects.
    */
   void setDefaultDepth(const int depth);

   /*!
    * @brief
    * Calculate the amount of memory needed to store the outeredge data
    * object, including object data and dynamically allocated data.
    */
   virtual size_t getSizeOfMemory(const hier::Box<DIM>& box) const;

   /*!
    * Return a boolean true value indicating that fine data for the outeredge quantity will
    * take precedence on coarse-fine interfaces.  See the OuteredgeVariable<DIM> class
    * header file for more information.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * Return true since the outeredge data index space extends beyond the interior of
    * patches.  That is, outeredge data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

   /*!
    * Return whether it is valid to copy this OuteredgeDataFactory to the 
    * supplied destination patch data factory.  It will return true if 
    * dst_pdf is EdgeDataFactory or OuteredgeDataFactory, false otherwise.
    */
   bool validCopyTo(
      const tbox::Pointer< hier::PatchDataFactory<DIM> >& dst_pdf) const;
   
private:
   int d_depth;

   hier::IntVector<DIM> *d_no_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuteredgeDataFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeDataFactory.C"
#endif

