//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/MemoryDatabase.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	An memomory database structure that stores (key,value) pairs in memory
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {



inline
MemoryDatabase::MemoryDatabase(const std::string& name)
{
   d_database_name = name;
}

inline
std::string
MemoryDatabase::getName() const
{
   return(d_database_name);
}

inline
bool MemoryDatabase::keyAccessed(const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return(keydata ? keydata->d_accessed : false);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
