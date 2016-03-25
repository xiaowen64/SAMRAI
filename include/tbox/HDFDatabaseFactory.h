//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2086 $
// Modified:	$LastChangedDate: 2008-03-28 15:10:07 -0700 (Fri, 28 Mar 2008) $
// Description:	An abstract base class for a HDFDatabaseFactory
//

#ifndef included_tbox_HDFDatabaseFactory
#define included_tbox_HDFDatabaseFactory

#include "SAMRAI_config.h"
#include "tbox/DatabaseFactory.h"

namespace SAMRAI {
   namespace tbox {

/**
 * @brief HDFDatabase factory.
 *
 * Builds a new HDFDatabase.
 */
class HDFDatabaseFactory : public DatabaseFactory
{	
   /**
    * Build a new Database object.
    */
   virtual Pointer<Database> allocate(const std::string& name);
};

}
}

#endif
