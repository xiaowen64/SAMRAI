//
// File:	RestartManager.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
void RestartManager::writeRestartFile(const string& root_dir_name)
{
   writeRestartFile(root_dir_name,0);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
