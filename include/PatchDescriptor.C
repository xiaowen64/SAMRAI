//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/PatchDescriptor.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Factory class for patch data objects that live on a patch
//

#ifndef included_hier_PatchDescriptor_C
#define included_hier_PatchDescriptor_C

#include "PatchDescriptor.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/Utilities.h"

#include <typeinfo>

#ifdef DEBUG_NO_INLINE
#include "PatchDescriptor.I"
#endif
namespace SAMRAI {
    namespace hier {

#define INDEX_UNDEFINED      (-1)

/*
*************************************************************************
*									*
* The constructor sets the max number of registered components to zero  *
* and allocates the factory and name arrays to the fixed length set     *
* by the SAMRAIManager utility.  The free list of indices               *
* is initialized to the full set of potentially used indices.           *
*									*
* The destructor clears the free index list and implicitly              *
* deallocates the arrays of name strings and factory pointers.          *  
*									*
*************************************************************************
*/

template<int DIM>  PatchDescriptor<DIM>::PatchDescriptor()
{
   const int max_num_patch_data_components_allowed = 
       tbox::SAMRAIManager::getMaxNumberPatchDataEntries(); 
   d_max_number_registered_components = 0;
   d_names.resizeArray(max_num_patch_data_components_allowed);
   d_factories.resizeArray(max_num_patch_data_components_allowed);
   for (int i = 0; i < max_num_patch_data_components_allowed; i++) {
      d_free_indices.appendItem(i);
   }
}

template<int DIM>  PatchDescriptor<DIM>::~PatchDescriptor()
{
   d_free_indices.clearItems();
}

/*
*************************************************************************
*									*
* Add the new factory to the list of patch data factories and assign    *
* it an integer index identifier.  Use a free list item if possible.    *
*									*
*************************************************************************
*/

template<int DIM> int PatchDescriptor<DIM>::definePatchDataComponent(
   const std::string& name,
   tbox::Pointer< PatchDataFactory<DIM> > factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !name.empty() );
   TBOX_ASSERT( !factory.isNull() );
#endif
   int ret_index = INDEX_UNDEFINED; 
   if (d_free_indices.isEmpty()) {
      TBOX_ERROR("PatchDescriptor::definePatchDataComponent error...\n"
                 << "No available patch data component indices left.\n"
                 << "Application must be restarted and size must be increased.\n" 
                 << "See tbox::SAMRAIManager utility for more information." << std::endl);
   } else {
      ret_index = d_free_indices.getFirstItem();
      d_free_indices.removeFirstItem();
      if (d_max_number_registered_components < ret_index+1) {
         d_max_number_registered_components = ret_index+1;
      }
      d_factories[ret_index] = factory;
      d_names    [ret_index] = name;
   }
   return(ret_index);
}

/*
*************************************************************************
*									*
* Remove the specified patch data factory index and place the index on	*
* the list of free indices.						*
*									*
*************************************************************************
*/

template<int DIM> void 
PatchDescriptor<DIM>::removePatchDataComponent(const int id)
{
   if ((id >= 0) && (id < d_max_number_registered_components)) {
      if (!d_names[id].empty()) {
         d_names[id] = std::string();
      }
      if (!d_factories[id].isNull()) {
         d_factories[id].setNull();
         d_free_indices.addItem(id);
      }
   }
}

/*
*************************************************************************
*									*
* Look up the factory by name; if no matching factory exists, then a	*
* pointer to null is returned.  The first matching factory is returned.	*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< PatchDataFactory<DIM> >
PatchDescriptor<DIM>::getPatchDataFactory(const std::string &name) const
{
   tbox::Pointer< PatchDataFactory<DIM> > factory = NULL;
   const int id = mapNameToIndex(name);
   if (id >= 0) {
      factory = d_factories[id];
   }
   return(factory);
}

/*
*************************************************************************
*									*
* Search the factory list for a match and return the associated		*
* factory.  If no match exists, return a negative identifier.		*
*									*
*************************************************************************
*/

template<int DIM> int 
PatchDescriptor<DIM>::mapNameToIndex(const std::string &name) const
{
   int ret_index = INDEX_UNDEFINED;
   int id = 0;
   while ( (ret_index == INDEX_UNDEFINED) && 
            (id < d_max_number_registered_components) ) {
      if (name == d_names[id]) {
         ret_index = id;
      }
      id++;
   }
   return(ret_index);
}

/*
*************************************************************************
*									*
* Print index, name, and factory data for the patch descriptor.		*
*									*
*************************************************************************
*/

template<int DIM> void PatchDescriptor<DIM>::printClassData(std::ostream& stream) const
{
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
   stream << "Printing PatchDescriptor<DIM> state ..." << std::endl;
   stream << "this = " << (PatchDescriptor<DIM>*)this << std::endl;
   stream << "d_max_number_registered_components = " 
          << d_max_number_registered_components << std::endl;
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
   for (int i = 0; i < d_max_number_registered_components; i++) {
      stream << "Patch Data Index=" << i << std::endl;
      if (!d_factories[i].isNull()) {
         stream << "   Patch Data Factory Name = " 
                << d_names[i] << std::endl;
         stream << "   Patch Data Factory = " 
                << typeid(*d_factories[i]).name() << std::endl;
      } else {
         stream << "   Patch Data Factory = NULL" << std::endl;
      }
   }
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

/*
*************************************************************************
*									*
* Return the maximum ghost cell width of all factories.          	*
*									*
*************************************************************************
*/

template<int DIM> IntVector<DIM>
PatchDescriptor<DIM>::getMaxGhostWidth() const
{
   IntVector<DIM> max_gcw(0);
   for (int i = 0; i < d_max_number_registered_components; i++) {
      if ( !d_factories[i].isNull() ) {
         max_gcw.max(d_factories[i]->getGhostCellWidth());
      }
   }
   return (max_gcw);
}

}
}
#endif
