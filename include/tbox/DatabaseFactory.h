//
// File:	$URL$
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2086 $
// Modified:	$LastChangedDate: 2008-03-28 15:10:07 -0700 (Fri, 28 Mar 2008) $
// Description:	An abstract base class for a DatabaseFactory
//

#ifndef included_tbox_DatabaseFactory
#define included_tbox_DatabaseFactory

#include "SAMRAI_config.h"

#include "tbox/Database.h"

namespace SAMRAI {
   namespace tbox {

/**
 * @brief Abstract base class factory used to build Database objects.
 * 
 * Used to build database objects.  For example, RestartManager 
 * may use a DatabaseFactory to build databases when creating
 * a restart file.
 */
class DatabaseFactory : public DescribedClass
{	
  public:
   /*
    * Build a new Database instance.
    */
   virtual Pointer<Database> allocate(const std::string& name) = 0;
};

}
}

#endif
