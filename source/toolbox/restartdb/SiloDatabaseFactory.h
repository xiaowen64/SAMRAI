//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building SiloDatabases
//

#ifndef included_tbox_SiloDatabaseFactory
#define included_tbox_SiloDatabaseFactory

#include "SAMRAI_config.h"
#include "tbox/DatabaseFactory.h"

namespace SAMRAI {
   namespace tbox {

/**
 * @brief SiloDatabase factory.
 *
 * Builds a new SiloDatabase.
 */
class SiloDatabaseFactory : public DatabaseFactory
{	
   /**
    * Build a new Database object.
    */
   virtual Pointer<Database> allocate(const std::string& name);
};

}
}

#endif
