//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/toolbox/inputdb/InputDatabase.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	An input database structure that stores (key,value) pairs
//

#ifndef included_tbox_InputDatabase
#define included_tbox_InputDatabase

#include "SAMRAI_config.h"

#include "tbox/MemoryDatabase.h"

namespace SAMRAI {
   namespace tbox {

/**
 * @brief Class InputDatabase stores (key,value) pairs in a hierarchical
 * database. 
 *
 * This is just another name for the MemoryDatabase. @see tbox::MemoryDatabase
 *
 * It is normally filled with data using a tbox::Parser (@see
 * tbox::Parser) and used to pass user supplied input from input files
 * to constructors for problem setup.
 *
 */

class InputDatabase : public MemoryDatabase
{
public:
   /**
    * The input database constructor creates an empty database with the
    * specified name.
    */
   InputDatabase(const std::string& name);

   /**
    * The input database destructor deallocates the data in the database.
    */
   virtual ~InputDatabase();
};

}
}

#endif
