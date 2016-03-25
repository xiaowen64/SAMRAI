//
// File:	InputDatabase.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 431 $
// Modified:	$Date: 2005-06-10 16:51:19 -0700 (Fri, 10 Jun 2005) $
// Description:	An input database structure that stores (key,value) pairs
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {



inline
InputDatabase::InputDatabase(const string& name)
{
   d_database_name = name;
}

inline
string
InputDatabase::getName() const
{
   return(d_database_name);
}

inline
bool InputDatabase::keyAccessed(const string& key)
{
   KeyData* keydata = findKeyData(key);
   return(keydata ? keydata->d_accessed : false);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
