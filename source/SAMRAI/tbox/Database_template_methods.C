/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for the SAMRAI database objects
 *
 ************************************************************************/

#ifndef included_tbox_Database_template_methods
#define included_tbox_Database_template_methods

#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

template<class TYPE>
void Database::putVector(
   const std::string& key,
   const std::vector<TYPE>& vector)
{
   unsigned int size = static_cast<int>(vector.size());
   putInteger(key + "_size", size);
   for (unsigned int i = 0; i < size; ++i) {
      const std::string index_str = tbox::Utilities::intToString(i);
      vector[i].putUnregisteredToDatabase(*this, key + "_" + index_str);
   }
}

template<class TYPE>
void Database::getVector(
   const std::string& key,
   std::vector<TYPE>& vector)
{
   size_t size = getInteger(key + "_size");
   for (unsigned int i = 0; i < size; ++i) {
      const std::string index_str = tbox::Utilities::intToString(i);
      vector[i].getFromDatabase(*this, key + "_" + index_str);
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
