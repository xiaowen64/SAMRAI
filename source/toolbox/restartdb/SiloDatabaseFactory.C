//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2196 $
// Modified:	$LastChangedDate: 2008-05-14 14:25:02 -0700 (Wed, 14 May 2008) $
// Description:	A factory for building SiloDatabases
//

#include "tbox/SiloDatabaseFactory.h"
#include "tbox/SiloDatabase.h"

namespace SAMRAI {
   namespace tbox {

#ifdef HAVE_SILO
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
#endif

}
}
