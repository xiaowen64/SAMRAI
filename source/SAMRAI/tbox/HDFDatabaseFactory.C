/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for a HDFDatabaseFactory 
 *
 ************************************************************************/

#include "SAMRAI/tbox/HDFDatabaseFactory.h"
#include "SAMRAI/tbox/HDFDatabase.h"

namespace SAMRAI {
namespace tbox {

/**
 * Build a new Database object.
 */
Pointer<Database> HDFDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_HDF5
   Pointer<HDFDatabase> database(new HDFDatabase(name));
   return database;

#else
   (void) name;
   TBOX_WARNING("HDF5DatabaseFactory: Cannot allocate an HDFDatabase.\n"
      << "SAMRAI was not configured with HDF.");
   return Pointer<Database>(NULL);
#endif
}

}
}
