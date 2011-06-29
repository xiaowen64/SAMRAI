/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A factory for building MemoryDatabases 
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabaseFactory.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

namespace SAMRAI {
namespace tbox {

/**
 * Build a new MemoryDatabase object.
 */
Pointer<Database> MemoryDatabaseFactory::allocate(
   const std::string& name) {
   Pointer<MemoryDatabase> database(new MemoryDatabase(name));
   return database;
}

}
}
