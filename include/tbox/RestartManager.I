//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/restartdb/RestartManager.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1769 $
// Modified:	$LastChangedDate: 2007-12-11 17:36:06 -0800 (Tue, 11 Dec 2007) $
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
