//
// File:	BasePatchHierarchy.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2003 The Regents of the University of California
// Revision:	$Revision: 551 $
// Modified:	$Date: 2005-08-17 11:15:27 -0700 (Wed, 17 Aug 2005) $
// Description:	An abstract base class of hierarchies
//

#ifndef included_hier_BasePatchHierarchy
#define included_hier_BasePatchHierarchy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BasePatchLevel
#include "BasePatchLevel.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif

namespace SAMRAI {
   namespace hier {


/*!
 * Class BasePatchHierarchy in a virtual base class that provides
 * an abstract interface for a patch hierarchy.  This allows higher-level
 * classes in SAMRAI and in applications that use SAMRAI to interface
 * with a hierarchy in an abstract manner without knowing whether it is a
 * PatchHierarchy representing a rectangular domain or if it is a
 * hierarchy that represents, for example, part of a multiblock domain.
 *
 * @see hier::PatchHierarchy
 */

template<int DIM> class BasePatchHierarchy
: public virtual tbox::DescribedClass,
  public tbox::Serializable
{
public:

   /*!
    * Default constructor
    */
   BasePatchHierarchy();

   /*!
    * Destructor for base patch hierarchy objects.
    */
   virtual ~BasePatchHierarchy();

   /*!
    * Return a pointer to the specified patch level.
    */
   virtual 
   tbox::Pointer<hier::BasePatchLevel<DIM> > getPatchLevel(const int l) const = 0;

   /*!
    * Returns true if the array of patch levels contains a patch level
    * finer than the specified patch level. Otherwise, false is returned.
    */
   virtual bool finerLevelExists(const int l) const = 0;

   /*!
    * Return the level number of the finest resolution patch level residing
    * in the hierarchy.
    */
   virtual int getFinestLevelNumber() const = 0;

   /*!
    * Return the number of levels that currently exist in the hierarchy.
    */
   virtual int getNumberLevels() const = 0;

   /**
    * Read in the entire hierarchy from the restart file.
    */
   virtual void getFromRestart(const int max_levels) = 0;

   /*!
    * Writes the state of the BasePatchHierarchy object and the PatchLevels
    * it contains to the database.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> database) = 0;

   /*!
    * Writes the state of the BasePatchHierarchy object and the PatchLevels
    * it contains to the database.  Only those patchdata corresponding to
    * the set bits in the ComponentSelector are written to the
    * specified database.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> database,
                              const ComponentSelector& patchdata_write_table) = 0;

private:


};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BasePatchHierarchy.C"
#endif

#endif
