//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/toolbox/inputdb/InputDatabase.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	An input database structure that stores (key,value) pairs
//

#include "tbox/InputDatabase.h"

namespace SAMRAI {
   namespace tbox {


/*
 * Constructor 
 */
InputDatabase::InputDatabase(const std::string& name) : MemoryDatabase(name)
{
}

/*
 * The virtual destructor deallocates database data.			
 */
InputDatabase::~InputDatabase()
{
}

}
}
