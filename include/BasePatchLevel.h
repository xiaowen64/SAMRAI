//
// File:	BasePatchLevel.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2003 The Regents of the University of California
// Revision:	$Revision: 551 $
// Modified:	$Date: 2005-08-17 11:15:27 -0700 (Wed, 17 Aug 2005) $
// Description:	An abstract base class for a level of the AMR hierarchy
//

#ifndef included_hier_BasePatchLevel
#define included_hier_BasePatchLevel

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
   namespace hier {


/**
 * Class BasePatchLevel in a virtual base class that provides
 * an abstract interface for a patch level.  This allows higher-level
 * classes in SAMRAI and in applications that use SAMRAI to interface
 * with a level in an abstract manner without knowing whether it is a
 * PatchLevel representing a level on a rectangular domain or if it is a
 * level that exists, for example, in a multiblock domain.
 *
 * @see hier::BasePatchHierarchy
 * @see hier::PatchLevel
 * @see hier::Patch
 */

template<int DIM> class BasePatchLevel
 : public tbox::DescribedClass
{
public:

   /**
    * Default constructor.  BasePatchLevel must be initialized before it can
    * be used.
    */
   BasePatchLevel();

   /**
    * The virtual destructor for patch level deallocates all patches.
    */
   virtual ~BasePatchLevel();

   /**
    * Return the number of this level in a hierarchy, or the number of
    * a hierarchy level matching the index space of this level. If this
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.  When the level is in a hierarchy, the return
    * value os the number of the level in the hierarchy.  See member 
    * function inHierarchy() below.
    */
   virtual int getLevelNumber() const = 0;

   /**
    * Allocate the specified component on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    */
   virtual
   void allocatePatchData(
      const int id,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL) = 0;
 
   /**
    * Allocate the specified components on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    */
   virtual
   void allocatePatchData(
      const hier::ComponentSelector& components,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL) = 0;
 
   /**
    * Deallocate the specified component on all patches.  This component 
    * will need to be reallocated before its next use.
    */
   virtual void deallocatePatchData(const int id) = 0;
 
   /**
    * Deallocate the specified components on all patches.  These components 
    * will need to be reallocated before their next use.
    */
   virtual void
   deallocatePatchData(const hier::ComponentSelector& components) = 0; 

   /**
    * Set the simulation time for the specified patch component.
    */
   virtual void setTime(const double timestamp, const int id) = 0;

   /**
    * Set the simulation time for the specified patch components.
    */
   virtual void setTime(
      const double timestamp,
      const hier::ComponentSelector& components) = 0;

   /**
    * Set the simulation time for all allocated patch components.
    */
   virtual void setTime(const double timestamp) = 0;

   /**
    * Return a const reference to the vector ratio between the index
    * space of this patch level and that of a reference level in AMR
    * hierarchy (typically, level zero).  Specifically, this is the
    * ratio passed to the constructor.
    */
   virtual const hier::IntVector<DIM>& getRatio() const = 0;

private:

};

}
}

#endif
