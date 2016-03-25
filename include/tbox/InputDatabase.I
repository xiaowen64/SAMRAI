//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/InputDatabase.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	An input database structure that stores (key,value) pairs
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {



inline
InputDatabase::InputDatabase(const std::string& name)
{
   d_database_name = name;
}

inline
std::string
InputDatabase::getName() const
{
   return(d_database_name);
}

inline
bool InputDatabase::keyAccessed(const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return(keydata ? keydata->d_accessed : false);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
