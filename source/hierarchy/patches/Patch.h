//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/Patch.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Patch container class for patch data objects
//

#ifndef included_hier_Patch
#define included_hier_Patch

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_stddef
#define included_stddef
#include <stddef.h>
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_hier_PatchGeometry
#include "PatchGeometry.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace hier {

/**
 * Class Patch<DIM> is a container for patch data objects defined over
 * some box region.  Each of the patch data objects that exist on the patch
 * are subclasses of the patch data pure virtual base class.  The construction
 * of the patch data objects on a patch are managed via the factory classes
 * maintained by the patch descriptor.
 *
 * By default, the patch constructor does not allocate space for the various
 * patch data components that live on the patch.  Individual components or sets
 * of components can be instantiated or destroyed via patch member functions.
 *
 * @see hier::Box
 * @see hier::PatchDescriptor
 * @see hier::PatchData
 * @see hier::PatchDataFactory
 * @see hier::PatchGeometry
 */

template<int DIM> class Patch  : public tbox::DescribedClass
{
public:
   /**
    * Allocate a patch container over the box, but do not instantiate 
    * any patch data components.
    */
   Patch(const Box<DIM>& box,
               tbox::Pointer< PatchDescriptor<DIM> > descriptor);

   /**
    * Virtual destructor for patch objects.
    */
   virtual ~Patch<DIM>();

   /**
    * Return the box over which this patch is defined.  The components on a
    * match are free to interpret this box as appropriate to the geometry of
    * the component (e.g., face, cell, or node centered).
    */
   const Box<DIM>& getBox() const;

   /**
    * Return the number of this patch.  In general, this number represents
    * the number of the patch in some patch level.
    */
   int getPatchNumber() const;

   /**
    * Set the number associated with this patch in the patch level.
    */
   void setPatchNumber(const int p);

   /**
    * Return the patch descriptor for this patch.  The patch descriptor
    * describes the patch data types that may exist on the patch.
    */
   tbox::Pointer< PatchDescriptor<DIM> > getPatchDescriptor() const;

   /**
    * Return a pointer to the patch data object associated with the specified
    * identifier.  Typically, this function should only be called to access
    * patch data objects that have already been explicitly allocated through 
    * a call to one of the patch allocation routines.  A NULL pointer will 
    * be returned if the patch data object is not allocated.  
    */
   tbox::Pointer< PatchData<DIM> > getPatchData(const int id) const;

   /**
    * Return a pointer to the patch data object associated with the specified
    * variable and context.  Typically, this function should only be called 
    * to access patch data objects that have already been explicitly allocated 
    * through a call to one of the patch allocation routines.  A NULL pointer 
    * will be returned if the patch data object is not allocated.
    */
   tbox::Pointer< PatchData<DIM> > getPatchData(
      const tbox::Pointer< Variable<DIM> > variable,
      const tbox::Pointer<VariableContext> context) const;

   /**
    * Set the patch data pointer associated with the specified identifier.
    * Note that this member function must be used with caution.  It can only
    * be called for patch data indices that were previously allocated through
    * one of the patch routines.  This member function does not check to see
    * whether the patch data types are already allocated, consistent, or
    * whether they have the same factory.  So, for example, a face centered
    * data type could be assigned to a location reserved for a cell centered
    * data type (with potentially different box and ghost cell width).
    */
   void setPatchData(const int id, tbox::Pointer< PatchData<DIM> > data);

   /**
    * Check whether the specified component has been allocated.
    */
   bool checkAllocated(const int id) const;

   /**
    * Return the amount of memory needed to allocate the specified component.
    */
   size_t getSizeOfPatchData(const int id) const;

   /**
    * Return the amount of memory needed to allocate the specified components.
    */
   size_t getSizeOfPatchData(const ComponentSelector& components) const;

   /**
    * Allocate the specified component on the patch.  If no memory
    * arena is specified, then the standard memory arena will be used.
    * Also, an assertion will result if the component is already allocated.
    * This provides a key bit of debugging information that may be useful
    * to application developers.
    */
   void allocatePatchData(const int id,
                          const double time = 0.0,
                          tbox::Pointer<tbox::Arena> pool = NULL);

   /**
    * Allocate the specified components on the patch.  If no memory
    * arena is specified, then the standard memory arena will be used.
    * Also, an assertion will result if any requested component is 
    * already allocated.  This provides a key bit of debugging information
    * that may be useful to application developers.
    */
   void allocatePatchData(const ComponentSelector& components,
                          const double time = 0.0,
                          tbox::Pointer<tbox::Arena> pool = NULL);

   /**
    * Deallocate the specified component.  This component will need to be
    * reallocated before its next use.
    */
   void deallocatePatchData(const int id);

   /**
    * Deallocate the specified components.  These components will need to be
    * reallocated before their next use.
    */
   void deallocatePatchData(const ComponentSelector& components);

   /**
    * Set the geometry specification for the patch domain.  This includes
    * patch boundary information and grid data.
    */
   void setPatchGeometry(tbox::Pointer< PatchGeometry<DIM> > geometry);

   /**
    * Return pointer to patch geometry object.
    */
   tbox::Pointer< PatchGeometry<DIM> > getPatchGeometry() const;

   /**
    * Set the simulation time for the specified patch component.
    */
   void setTime(const double timestamp, const int id);

   /**
    * Set the simulation time for the specified patch components.
    */
   void setTime(const double timestamp,
                const ComponentSelector& components);

   /**
    * Set the simulation time for all allocated patch components.
    */
   void setTime(const double timestamp);

   /**
    * Get the level number of the patch level where this patch resides.
    * Note that this value can be a valid level number (i.e., >=0) even
    * when the level is not in a hierarchy.  In this case the level
    * number will be the number of the hierarchy level matching the
    * index space of the patch level holding this patch.  If the patch
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.  See member function inHierarchy() below.
    */
   int getPatchLevelNumber() const;

   /**
    * Set the patch level number to the number of a level in the hierarchy
    * when the index space of the level owning this patch aligns with the
    * index space of some valid hierarchy level.   The default level
    * number is -1.
    */
   void setPatchLevelNumber(const int level_number);

   /**
    * Return true when the level holding this patch resides in a hierarchy
    * and false otherwise.
    */
   bool inHierarchy() const;

   /**
    * Set to true if the level holding this patch resides in a hierarchy;
    * false otherwise.  The default setting is false.
    */
   void setPatchInHierarchy(bool in_hierarchy);

   /**
    * Checks that class version and restart file version are equal.  
    * If so, read in the state of the patch and create all patch data
    * objects indicated in the component selector from the database
    *
    * Assertions: check that database is a non-null Pointer,
    * that data retrieved from the database are of the type
    * expected, and that the patch_number read in from the database
    * matches the patch number assigned to this Patch.
    *
    * Important Note: a warning will be printed to the log file if 
    * some patch data components that were requested through the 
    * component_selector are not found in the database.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database> database,
                        ComponentSelector component_selector);

   /**
    * Write out the class version number and the state of the patch
    * object.  Finally, patch data objects specified in the component
    * selector are written to the specified database.
    * 
    * Assertions: check that database is a non-null Pointer.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      const ComponentSelector& patchdata_write_table);


   /*!
    * @brief Print a patch (for debugging).
    *
    * Depth is kept for consistency with other recursivePrint methods,
    * but is not used because this is the lowest level of recursion
    * currently supported.
    */
   int recursivePrint( std::ostream &os ,
                       const std::string &border=std::string() ,
                       unsigned short depth=0 ) const;

   /**
    * Output patch information (box and number of components).
    */
   template<int D> 
      friend std::ostream& operator<<(std::ostream& s, const Patch<D>& patch);

private:
   Patch(const Patch<DIM>&);	// not implemented
   void operator=(const Patch<DIM>&);	// not implemented

   Box<DIM>                                   d_box;
   tbox::Pointer< PatchDescriptor<DIM> >         d_descriptor;
   tbox::Pointer< PatchGeometry<DIM> >           d_patch_geometry;
   tbox::Array< tbox::Pointer< PatchData<DIM> > > d_patch_data;
   int                                         d_patch_number;
   int                                         d_patch_level_number;
   bool                                        d_patch_in_hierarchy;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "Patch.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "Patch.C"
#endif
