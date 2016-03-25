/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for patch data objects that live on a patch
 *
 ************************************************************************/

#ifndef included_hier_PatchDescriptor
#define included_hier_PatchDescriptor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/List.h"

#include <string>
#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class PatchDescriptor maintains a collection of patch data
 * factories and associated names that describes how patch data entries are
 * constructed on each patch in an AMR hierarchy.  The factory mechanism is
 * used to create new instances of concrete patch data objects without knowing
 * their actual types.  See the Design Patterns book by Gamma {\em et al.}
 * for more details about the Abstract Factory pattern.  Generally, a PatchDescriptor
 * object is intended to be shared among all patches (which are distributed across
 * processors) so that they store patch data objects in the same way.
 *
 * Patch data factory objects (and associated names) are explicitly added to the
 * PatchDescriptor using the definePatchDataComponent() member function.  This function
 * returns an integer index that can be used to identify the corresponding patch data
 * on a a patch.  Factories can be removed from the PatchDescriptor using the
 * removePatchDataComponent() member function, which returns the integer index associated
 * with the removed factory to a "free list" so that it can be used again.  At any time,
 * the valid range of indices is >= 0 and < getMaxNumberRegisteredComponents().
 *
 * Note that the SAMRAIManager utility establishes a maximum number of patch data
 * object that may live on a Patch object which, for consistency, must be the same as
 * the number of patch data factories a PatchDescriptor will hold.   See the documentation
 * of the SAMRAIManager utility for information about changing this maximum value.
 *
 * @see tbox::SAMRAIManager
 * @see hier::PatchDataFactory
 * @see hier::PatchDataData
 * @see hier::Patch
 */

class PatchDescriptor:public tbox::DescribedClass
{
public:
   /*!
    * Constructor for a patch descriptor initializes the
    * descriptor to hold zero patch data factory entries.
    */
   PatchDescriptor();

   /*!
    * The virtual destructor for a patch descriptor deallocates the
    * internal data structures.
    */
   virtual ~PatchDescriptor();

   /*!
    * Add a new patch data factory and name std::string identifier to the patch
    * descriptor.  The factory will be given the specified name which must be
    * unique for the mapNameToIndex() function to execute as expected. However,
    * there is no internal checking done to ensure that names are unique.
    *
    * @return int index assigned to given patch data factory in patch descriptor.
    *
    * @param name      std::string name to be associated in name list with given factory,
    *                  which must be non-empty when assertion checking is active.
    * @param factory   pointer to factory to add to patch descriptor, which must
    *                  be non-null when assertion checking is active.
    */
   int
   definePatchDataComponent(
      const std::string& name,
      tbox::Pointer<PatchDataFactory> factory);

   /*!
    * Deallocate the patch data factory in the patch descriptor identified by the
    * given index.  The index may be assigned to another factory in the future.
    * However, index will be invalid as a patch data index until it is re-allocated
    * by the definePatchDataComponent() member function.  An invalid id value
    * passed to this function is silently ignored.
    *
    * @param id      int index of factory to remove from patch descriptor.
    */
   void
   removePatchDataComponent(
      int id);

   /*!
    * Retrieve a patch data factory by integer index identifier.  The identifier
    * is the one previously returned by definePatchDataComponent().  Note that the
    * factory pointer will be null if the index is is not currently assigned.
    *
    * @return pointer to patch data factory assigned to given index.
    *
    * @param id      int index of factory to return, which must be >= 0 and
    *                < the return value of getMaxNumberRegisteredComponents();
    */
   tbox::Pointer<PatchDataFactory>
   getPatchDataFactory(
      int id) const;

   /*!
    * Retrieve a patch data factory by name std::string identifier.  Recall that
    * uniqueness of names is not strictly enforced. So if more than one
    * factory matches the given name, then only one of them is returned.  If no
    * matching factory is found, then a null pointer is returned.
    *
    * @return pointer to patch data factory assigned to given name.
    *
    * @param name    std::string name of factory.
    */
   tbox::Pointer<PatchDataFactory>
   getPatchDataFactory(
      const std::string& name) const;

   /*!
    * Get the maximum number of components currently known to the patch
    * descriptor.  That is, this number indicates the largest number of
    * components that have been registered with the descriptor via the
    * definePatchDataComponent() function, which is equal to the largest
    * known patch data component index + 1.  Note that the total number of
    * registered components is reduced by calls to removePatchDataComponent(),
    * but the max number remains the same when components are removed.
    * In that case, the corresponding indices are placed on a list of "free"
    * values to be re-used in subsequent calls to definePatchDataComponent().
    *
    * @return largest index assigned to this point.
    */
   int
   getMaxNumberRegisteredComponents() const;

   /*!
    * Lookup a factory by std::string name and return its integer index identifier.
    * Note that more than one factory may have the same name.  In this case, the
    * identifier of one of the factories is chosen.  If no matching factory is found,
    * then an invalid negative index is returned.
    */
   int
   mapNameToIndex(
      const std::string& name) const;

   /*!
    * Lookup a factory by identifier and return its name.
    */
   const std::string&
   mapIndexToName(
      const int id) const;

   /*!
    * Return the IntVector indicating the maximum ghost cell width of all registered
    * patch data components for the provided dimension.
    *
    * If no components have been registered returns the value set by
    * setMinGhostWidth(), which is zero by default.
    *
    * @param dim Dimension
    */
   IntVector
   getMaxGhostWidth(
      const tbox::Dimension& dim) const;

   /*!
    * @brief Set a minimum value on the value returned by
    * getMaxGhostWidth().
    *
    * This method allows users to specify a mininum value returned by
    * getMaxGhostWidth().  The default minimum is zero.  This value
    * can be used as a substitute for data that is not yet registered
    * with the PatchDescriptor and therefore cannot be reflected in
    * getMaxGhostWidth().
    *
    * The dimension associated with the set value is taken to be @c
    * min_value.getDim().
    */
   void
   setMinGhostWidth(
      const IntVector& min_value);

   /*!
    * Print patch descriptor data to given output stream (plog by default).
    */
   virtual void
   printClassData(
      std::ostream& stream = tbox::plog) const;

private:
   /*
    * Static integer constant describing value of an undefined index.
    */
   static const int INDEX_UNDEFINED;

   PatchDescriptor(
      const PatchDescriptor&);                  // not implemented
   void
   operator = (
      const PatchDescriptor&);                  // not implemented

   int d_max_number_registered_components;
   tbox::Array<std::string> d_names;
   tbox::Array<tbox::Pointer<PatchDataFactory> > d_factories;
   tbox::List<int> d_free_indices;

   /*!
    * @brief Value set by setMinGhostWidth().
    */
   std::vector<IntVector> d_min_gcw;

};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/PatchDescriptor.I"
#endif
#endif
