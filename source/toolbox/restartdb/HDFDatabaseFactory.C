//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/toolbox/database/Database.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2086 $
// Modified:	$LastChangedDate: 2008-03-28 15:10:07 -0700 (Fri, 28 Mar 2008) $
// Description:	An abstract base class for a HDFDatabaseFactory
//

#include "tbox/HDFDatabaseFactory.h"
#include "tbox/HDFDatabase.h"

namespace SAMRAI {
   namespace tbox {

/**
 * Build a new Database object.
 */
Pointer<Database> HDFDatabaseFactory::allocate(const std::string& name) {
   Pointer<HDFDatabase> database = new HDFDatabase(name);
   return database;
}

}
}
