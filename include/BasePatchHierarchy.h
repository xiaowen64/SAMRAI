//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/BasePatchHierarchy.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	An abstract base class of hierarchies
//

#ifndef included_hier_BasePatchHierarchy
#define included_hier_BasePatchHierarchy

#include "SAMRAI_config.h"
#include "BasePatchLevel.h"
#include "IntVector.h"
#include "PatchDescriptor.h"
#include "ProcessorMapping.h"
#include "tbox/DescribedClass.h"
#include "tbox/Serializable.h"

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
   virtual int getNumberOfLevels() const = 0;

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
