/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A factory for building SiloDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/SiloDatabaseFactory.h"
#include "SAMRAI/tbox/SiloDatabase.h"

namespace SAMRAI {
namespace tbox {

#ifdef HAVE_SILO
/**
 * Build a new SiloDatabase object.
 */
Pointer<Database> SiloDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_SILO
   Pointer<SiloDatabase> database(new SiloDatabase(name));
   return database;

#else
   return NULL;

#endif
}
#endif

}
}
