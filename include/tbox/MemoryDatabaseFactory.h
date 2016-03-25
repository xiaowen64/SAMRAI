//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building MemoryDatabases
//

#ifndef included_tbox_MemoryDatabaseFactory
#define included_tbox_MemoryDatabaseFactory

#include "SAMRAI_config.h"
#include "tbox/DatabaseFactory.h"

namespace SAMRAI {
   namespace tbox {

/**
 * @brief MemoryDatabase factory.
 *
 * Builds a new MemoryDatabase.
 */
class MemoryDatabaseFactory : public DatabaseFactory
{	
   /**
    * Build a new Database object.
    */
   virtual Pointer<Database> allocate(const std::string& name);
};

}
}

#endif
