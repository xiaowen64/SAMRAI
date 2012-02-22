/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An input manager singleton class that parses input files
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Return whether or not the manager contains an valid input database.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool InputManager::inputDatabaseExists()
{
   return s_input_db;
}

/*
 *************************************************************************
 *
 * Accessor method for InputManger's root input database.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
boost::shared_ptr<Database> InputManager::getInputDatabase()
{
   return s_input_db;
}

}
}
