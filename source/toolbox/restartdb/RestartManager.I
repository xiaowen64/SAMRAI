//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/restartdb/RestartManager.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2169 $
// Modified:	$LastChangedDate: 2008-04-30 16:24:39 -0700 (Wed, 30 Apr 2008) $
// Description:	An restart manager singleton class 
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* Destructor								*
*									*
*************************************************************************
*/
inline
RestartManager::~RestartManager()
{
}

/*
*************************************************************************
*									*
* Returns d_is_from_restart.						*
*									*
*************************************************************************
*/

inline
bool RestartManager::isFromRestart()
{
   return d_is_from_restart;
}

/*
*************************************************************************
*									*
* Returns the root of the restart database which was opened at the 	*
* start of the simulation.						*
*									*
*************************************************************************
*/
inline
Pointer<Database> RestartManager::getRootDatabase()
{
   return d_database_root;
}

/*
*************************************************************************
*									*
* Sets the root database                                                *
*									*
*************************************************************************
*/
inline
void RestartManager::setRootDatabase(Pointer<Database> database)
{
   if(database.isNull()) {
      d_database_root = NULL;
      d_is_from_restart = false; 
   } else {
      d_database_root = database;
      d_is_from_restart = true; 
   }
}

/*
*************************************************************************
*									*
* Sets the database factory                                             *
*									*
*************************************************************************
*/
inline
void RestartManager::setDatabaseFactory(Pointer<DatabaseFactory> database_factory)
{
   d_database_factory = database_factory;
}


/*
*************************************************************************
*									*
* writeRestartFile							*
*									*
*************************************************************************
*/
inline
void RestartManager::writeRestartFile(const std::string& root_dir_name)
{
   writeRestartFile(root_dir_name,0);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
