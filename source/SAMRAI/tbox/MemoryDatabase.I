/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An memomory database structure that stores (key,value) pairs in memory 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
MemoryDatabase::MemoryDatabase(
   const std::string& name)
{
   d_database_name = name;
}

SAMRAI_INLINE_KEYWORD
std::string
MemoryDatabase::getName() const
{
   return d_database_name;
}

SAMRAI_INLINE_KEYWORD
bool MemoryDatabase::keyAccessed(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_accessed : false;
}

}
}
