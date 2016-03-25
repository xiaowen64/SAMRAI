//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/PatchData.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Abstract base class for patch data objects
//

#ifndef included_hier_PatchData_C
#define included_hier_PatchData_C

#include "PatchData.h"


#include "tbox/Utilities.h"

#define HIER_PATCH_DATA_VERSION (2)

#ifdef DEBUG_NO_INLINE
#include "PatchData.I"
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchData<DIM>::PatchData(const Box<DIM>& domain,
                                 const IntVector<DIM>& ghosts)
{
   d_box       = domain;
   d_ghosts    = ghosts;
   d_ghost_box = Box<DIM>::grow(domain, ghosts);
   d_timestamp = 0.0;
}

template<int DIM>  PatchData<DIM>::~PatchData()
{
}

/*
*************************************************************************
*                                                                       *
* Checks that clas and restart file version number are same.  If so, 	*
* reads in data members common to all patch data and then invoke	*
* getSpecializedFromDatabase() to read in data particular to the	*
* specific derived class.						*
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchData<DIM>::getFromDatabase(tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("HIER_PATCH_DATA_VERSION");
   if (ver != HIER_PATCH_DATA_VERSION) {
      TBOX_ERROR("PatchData<DIM>::getFromDatabase() error...\n"
              << "  Restart file version different than class version" << std::endl);
   }

   d_box = database->getDatabaseBox("d_box");
   d_ghost_box = database->getDatabaseBox("d_ghost_box");
   database->getIntegerArray("d_ghosts", (int*)d_ghosts, DIM);
   d_timestamp = database->getDouble("d_timestamp");

   getSpecializedFromDatabase(database);
}

/*
*************************************************************************
*                                                                       *
* Write out data members common to all patch data and then invoke	*
* putSpecializedToDatabase() to write out data particular to the	*
* specific derived class.						*
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchData<DIM>::putToDatabase(tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif
   
   database->putInteger("HIER_PATCH_DATA_VERSION", HIER_PATCH_DATA_VERSION);
   database->putDatabaseBox("d_box", d_box);
   database->putDatabaseBox("d_ghost_box", d_ghost_box);
   database->putDouble("d_timestamp", d_timestamp);
   database->putIntegerArray("d_ghosts", (int*)d_ghosts, DIM);

   putSpecializedToDatabase(database);
}

}
}
#endif
