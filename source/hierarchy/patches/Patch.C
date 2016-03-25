//
// File:	Patch.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 407 $
// Modified:	$Date: 2005-06-01 09:49:54 -0700 (Wed, 01 Jun 2005) $
// Description:	Patch container class for patch data objects
//

#ifndef included_hier_Patch_C
#define included_hier_Patch_C

#include <typeinfo>
#include <string>
using namespace std;

#include "Patch.h"

#include "tbox/ArenaManager.h"
#define included_String
#include "tbox/Utilities.h"
#include "PatchDataFactory.h"
#include "VariableDatabase.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#define HIER_PATCH_VERSION (2)

#ifdef DEBUG_NO_INLINE
#include "Patch.I"
#endif
namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*									*
* Allocate a patch container but do not instantiate any components.	*
*									*
*************************************************************************
*/

template<int DIM>  Patch<DIM>::Patch(const Box<DIM>& box,
                         tbox::Pointer< PatchDescriptor<DIM> > descriptor)
{
   d_box          = box;
   d_descriptor   = descriptor;
   d_patch_number = -1;
   d_patch_level_number = -1;
   d_patch_in_hierarchy = false;
   d_touches_regular_boundary = false;
   d_touches_periodic_boundary = false;
}

/*
*************************************************************************
*									*
* The virtual destructor does nothing; all memory deallocation is	*
* managed automatically by the pointer and array classes.		*
*									*
*************************************************************************
*/

template<int DIM>  Patch<DIM>::~Patch()
{
}

/*
*************************************************************************
*									*
* Calculate the amount of memory space required to allocate the		*
* specified component(s).  This information can then be used by a	*
* fixed-size memory allocator.						*
*									*
*************************************************************************
*/

template<int DIM> size_t Patch<DIM>::getSizeOfPatchData(const int id) const
{
   return(d_descriptor->getPatchDataFactory(id)->getSizeOfMemory(d_box));
}

template<int DIM> size_t
Patch<DIM>::getSizeOfPatchData(const ComponentSelector& components) const
{
   size_t size = 0;
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();

   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i)) {
         size += d_descriptor->getPatchDataFactory(i)->getSizeOfMemory(d_box);
      }
   }

   return(size);
}

/*
*************************************************************************
*									*
* Allocate the specified patch data object(s) on the patch.  If no	*
* arena is specified, then the standard memory arena will be used.	*
*									*
*************************************************************************
*/

template<int DIM> void Patch<DIM>::allocatePatchData(const int id,
                                    const double time,
                                    tbox::Pointer<tbox::Arena> pool)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((id >= 0) && (id < ncomponents));
#endif

   if (ncomponents > d_patch_data.getSize()) {
      d_patch_data.resizeArray(ncomponents);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_patch_data[id].isNull());
#endif

   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   d_patch_data[id] =
      d_descriptor->getPatchDataFactory(id)->allocate(d_box, pool);
   d_patch_data[id]->setTime(time);
}
   
template<int DIM> void Patch<DIM>::allocatePatchData(const ComponentSelector& components,
                                    const double time,
                                    tbox::Pointer<tbox::Arena> pool)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();
   if (ncomponents > d_patch_data.getSize()) {
      d_patch_data.resizeArray(ncomponents);
   }

   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i)) {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_patch_data[i].isNull());
#endif
         d_patch_data[i] =
            d_descriptor->getPatchDataFactory(i)->allocate(d_box, pool);
         d_patch_data[i]->setTime(time);
      }
   }
}
 
/*
*************************************************************************
*									*
* Deallocate (or set to null) the specified component(s).		*
*									*
*************************************************************************
*/

template<int DIM> void Patch<DIM>::deallocatePatchData(const int id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((id >= 0) && (id < d_descriptor->getMaxNumberRegisteredComponents()));
#endif
   if (id < d_patch_data.getSize()) {
      d_patch_data[id].setNull();
   }
}

template<int DIM> void Patch<DIM>::deallocatePatchData(const ComponentSelector& components)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i)) {
         d_patch_data[i].setNull();
      }
   }
}

/*
*************************************************************************
*									*
* Set the time stamp for the specified components in the patch.		*
*									*
*************************************************************************
*/

template<int DIM> void Patch<DIM>::setTime(const double timestamp,
                          const ComponentSelector& components)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i) && !d_patch_data[i].isNull()) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

template<int DIM> void Patch<DIM>::setTime(const double timestamp)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (!d_patch_data[i].isNull()) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

/*
*************************************************************************
*									*
* Checks that class and restart file version numbers are equal.  If so, *
* reads in data from database and have each patch_data item read 	*
* itself in from the database						*
*									*
*************************************************************************
*/

template<int DIM> void Patch<DIM>::getFromDatabase(
   tbox::Pointer<tbox::Database> database,
   ComponentSelector component_selector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif

   int ver = database->getInteger("HIER_PATCH_VERSION");
   if (ver != HIER_PATCH_VERSION) {
      TBOX_ERROR("Patch<DIM>::getFromDatabase() error...\n"
         << "   Restart file version different than class version" << endl);
   }

   d_box = database->getDatabaseBox("d_box");
   int patch_number = database->getInteger("d_patch_number");

   if (patch_number != d_patch_number) {
      TBOX_ERROR("Patch<DIM>::getFromDatabase() error...\n"
              << "    patch number " << patch_number 
              << " read from database does not "
              << "match patch number of patch " << d_patch_number
              << endl);
   }

   d_patch_level_number = database->getInteger("d_patch_level_number");
   d_patch_in_hierarchy = database->getBool("d_patch_in_hierarchy");

   d_patch_data.resizeArray(d_descriptor->getMaxNumberRegisteredComponents());

   int namelist_count = database->getInteger("patch_data_namelist_count");
   tbox::Array<string> patch_data_namelist;
   if ( namelist_count ) {
      patch_data_namelist = database->getStringArray("patch_data_namelist");
   }

   for (int i = 0; i < patch_data_namelist.getSize(); i++) {
      string patch_data_name;
      int patch_data_index;
      tbox::Pointer<tbox::Database> patch_data_database;
      tbox::Pointer< PatchDataFactory<DIM> > patch_data_factory;

      patch_data_name  = patch_data_namelist[i];

      if (database->isDatabase(patch_data_name)) {
         patch_data_database = database->getDatabase(patch_data_name);
      } else {
         TBOX_ERROR("Patch<DIM>::getFromDatabase() error...\n"
                 << "   patch data" << patch_data_name 
                 << " not found in database" << endl);
      }

      patch_data_index = d_descriptor->
                         mapNameToIndex(patch_data_name);

      if ( (patch_data_index >= 0) && 
           (component_selector.isSet(patch_data_index)) ) {
         patch_data_factory = d_descriptor->
                                 getPatchDataFactory(patch_data_index);
         d_patch_data[patch_data_index] = 
                                 patch_data_factory->allocate(d_box);
         d_patch_data[patch_data_index]->getFromDatabase(patch_data_database);

         component_selector.clrFlag(patch_data_index);
      }
   }

   ComponentSelector all_unset(false);
   if (component_selector != all_unset) {
      TBOX_WARNING("Patch<DIM>::getFromDatabase() warning...\n"
                << "   Some requested patch data components not "
                << "found in database" << endl);
   }

}

/*
*************************************************************************
*									*
* Write out the class version number to database.  Then,		*
* writes out data to database and have each patch_data item write 	*
* itself out to the database.  The following data                       *
* members are written out: d_box, d_patch_number, d_patch_level_number, *
* d_patch_in_hierarchy, d_patch_data[].                                 *
* The database key for all data members is identical to the             *
* name of the data member except for the d_patch_data.  These have      *
* keys of the form "variable##context" which is the form that they      *
* are stored by the patch descriptor.  In addition a list of the        *
* patch_data names ("patch_data_namelist") and the number of patch data *
* items saved ("namelist_count") are also written to the database.      *
* The patchdata_write_table determines which patchdata are written to   *
* the database.                                                         *
*									*
*************************************************************************
*/
template<int DIM> void Patch<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector& patchdata_write_table)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif
   int i;

   database->putInteger("HIER_PATCH_VERSION", HIER_PATCH_VERSION);
   database->putDatabaseBox("d_box", d_box);
   database->putInteger("d_patch_number", d_patch_number);
   database->putInteger("d_patch_level_number", d_patch_level_number);
   database->putBool("d_patch_in_hierarchy", d_patch_in_hierarchy);

   int namelist_count = 0;
   for (i = 0; i < d_patch_data.getSize(); i++) {
      if ( patchdata_write_table.isSet(i) && checkAllocated(i) ) {
         namelist_count++;
      }
   }

   string patch_data_name;
   tbox::Pointer<tbox::Database> patch_data_database;
   tbox::Array<string> patch_data_namelist(namelist_count);
   namelist_count = 0;
   for (i = 0; i < d_patch_data.getSize(); i++) {
      if ( patchdata_write_table.isSet(i) && checkAllocated(i) ) {
         patch_data_namelist[namelist_count++] = 
            patch_data_name = d_descriptor->mapIndexToName(i);
         patch_data_database = database->putDatabase(patch_data_name);
         (d_patch_data[i])->putToDatabase(patch_data_database);
      }
   }

   database->putInteger("patch_data_namelist_count", namelist_count);
   if ( namelist_count > 0 ) {
      database->putStringArray("patch_data_namelist", patch_data_namelist);
   }
}

/*
*************************************************************************
*									*
* Print information about the patch.					*
*									*
*************************************************************************
*/


template<int DIM> int Patch<DIM>::recursivePrint( ostream &os ,
                                                  const string &border ,
                                                  unsigned short depth ) const
{
  os << border
     << d_box
     << "\tdims: " << d_box.numberCells(0)
     ;
  for ( int i=1; i<DIM; ++i ) {
     os << " X " << d_box.numberCells(i);
  }
  os << "\tsize: " << d_box.size()
     << "\n";
  return 0;
}

template<int DIM> ostream& operator<<(ostream& s, const Patch<DIM>& patch)
{
   s << "Patch<DIM>::box = " << patch.d_box << endl << flush;
   s << "Patch<DIM>::patch_number = " << patch.d_patch_number 
     << endl << flush;
   s << "Patch<DIM>::patch_level_number = " << patch.d_patch_level_number 
     << endl << flush;
   s << "Patch<DIM>::patch_in_hierarchy = " << patch.d_patch_in_hierarchy 
     << endl << flush;
   s << "Patch<DIM>::number_components = " << patch.d_patch_data.getSize()
     << endl << flush;
   const int ncomponents = patch.d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      s << "Component(" << i << ")=";
      if (patch.d_patch_data[i].isNull()) {
         s << "NULL\n";
      } else {
         s << typeid(*patch.d_patch_data[i]).name()
           << " [GCW=" << patch.d_patch_data[i]->getGhostCellWidth() << "]\n";
      }
   }
   return(s);
}

}
}
#endif
