//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building MemoryDatabases
//

#include "tbox/MemoryDatabaseFactory.h"
#include "tbox/MemoryDatabase.h"

namespace SAMRAI {
   namespace tbox {

/**
 * Build a new MemoryDatabase object.
 */
Pointer<Database> MemoryDatabaseFactory::allocate(const std::string& name) {
   Pointer<MemoryDatabase> database = new MemoryDatabase(name);
   return database;
}

}
}
