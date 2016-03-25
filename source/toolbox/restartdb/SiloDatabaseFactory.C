//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building SiloDatabases
//

#include "tbox/SiloDatabaseFactory.h"
#include "tbox/SiloDatabase.h"

namespace SAMRAI {
   namespace tbox {

/**
 * Build a new SiloDatabase object.
 */
Pointer<Database> SiloDatabaseFactory::allocate(const std::string& name) {
#ifdef HAVE_SILO
   Pointer<SiloDatabase> database = new SiloDatabase(name);
   return database;
#else
   return NULL;
#endif
}

}
}
