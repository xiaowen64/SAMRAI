/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An memomory database structure that stores (key,value) pairs
 *                in memory
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
bool MemoryDatabase::keyAccessed(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_accessed : false;
}

/*
 *************************************************************************
 *
 * Print unused database keys to the specified output stream.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void MemoryDatabase::printUnusedKeys(
   std::ostream& os) const
{
   printDatabase(os, 0, PRINT_UNUSED);
}

/*
 *************************************************************************
 *
 * Print default database keys to the specified output stream.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void MemoryDatabase::printDefaultKeys(
   std::ostream& os) const
{
   printDatabase(os, 0, PRINT_DEFAULT);
}

/*
 *************************************************************************
 *
 * Search the current database for a matching key.  If found, delete
 * that key and return true.  If the key does not exist, then return
 * false.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool MemoryDatabase::deleteKeyIfFound(
   const std::string& key)
{
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      if (i().d_key == key) {
         d_keyvalues.removeItem(i);
         return true;
      }
   }
   return false;
}

/*
 *************************************************************************
 *
 * Find the key data associated with the specified key and return a
 * pointer to the record.  If no such key data exists, then return NULL.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
MemoryDatabase::KeyData *
MemoryDatabase::findKeyData(
   const std::string& key)
{
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      if (key == i().d_key) return &i();
   }
   return NULL;
}

/*
 *************************************************************************
 *
 * Indent the output stream by the specified indentation factor.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void MemoryDatabase::indentStream(
   std::ostream& os,
   const long indent)
{
   for (int i = 0; i < indent; i++) {
      os << " ";
   }
}

}
}
